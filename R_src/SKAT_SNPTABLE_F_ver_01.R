# library path
library_load <- function(){
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(tidyverse);library(parallel)
  library(tidyselect);library(magrittr);library(SKAT);
  
  # tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/tool/:")
}

## SKAT-O FUNCTION
{
  # function
  geneset_load_SKAT <- function(){
    print("Geneset load_0428!")
    return_list <- list()
    geneset_merge <- fread(file = "/home/jinoo/skat-o/SKAT_data/1 gene sets 190428_exclude_ARX_for SKAT_O.txt", header = T, 
                           sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
    
    one_hot <- lapply(X=select(geneset_merge, -HGNC, -M2, -'M2-GBA-LRRK2', -'M2_Brain_gene', -'M2_NOT_brain'), 
                      FUN = function(x){ifelse(nchar(x) > 0, 1, 0)}) %>% bind_cols()
    geneset_onehot <- cbind(GENE=geneset_merge$HGNC, one_hot) %>% mutate(GENE = as.character(GENE)) %>% 
      cbind(M2 = as.numeric(geneset_merge$M2)) %>% cbind('M2-GBA-LRRK2'= as.numeric(geneset_merge$`M2-GBA-LRRK2`)) %>%
      cbind('M2_Brain_gene' = as.numeric(geneset_merge$`M2_Brain_gene`)) %>% 
      cbind('M2_NOT_brain'= as.numeric(geneset_merge$`M2_NOT_brain`))
    
    geneset_onehot <- geneset_onehot %>% mutate(., 'M2_Brain_gene-LRRK2-GBA' = M2_Brain_gene)
    geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[554]] <- 0
    geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[755]] <- 0
    
    for(geneset in 2:ncol(geneset_onehot)){
      for(k in 1:nrow(geneset_onehot)){
        if(geneset_onehot[k, geneset] >= 1){
          geneset_onehot[k,geneset] <- geneset_onehot[k,1]
        } else{
          geneset_onehot[k,geneset] <- NA
        }
      }
    }
    
    # CHD, DEG adding
    
    CHD <- fread(file = "/home/jinoo/skat-o/SKAT_data/CHD_0424.txt", header = T, 
                 sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
    CHD[CHD == "05-Mar"] <- "MARCH5" #### 0718 HGNC comfirm
    
    CHD <- CHD %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.)), CHD = NA))
    geneset_onehot <- bind_cols(geneset_onehot, CHD)
    
    DEG <- fread(file = "/home/jinoo/skat-o/SKAT_data/DEG_geneset.txt", header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    meta_DEG <- fread(file = "/home/jinoo/skat-o/SKAT_data/$Table 3 and 4. result_down_intersect_190520_0603.txt", header = T, data.table = F)
    
    
    DEG1 <- DEG$DEG1 %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG1) <- "DEG1"
    DEG2 <- DEG$DEG1YJK[DEG$DEG1YJK != ""] %>% as_tibble() %>% 
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG2) <- "DEG1YJK"
    DEG3 <- DEG$DEG1Mito[DEG$DEG1Mito != ""] %>% as_tibble() %>% 
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG3) <- "DEG1Mito"
    meta_DEG <- meta_DEG$SYMBOL[meta_DEG$SYMBOL != ""] %>% as_tibble() %>%
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(meta_DEG) <- "META_DEG"
    
    geneset_onehot <- geneset_onehot %>% bind_cols(., DEG1) %>% bind_cols(., DEG2) %>% bind_cols(., DEG3) %>%
      bind_cols(., meta_DEG)
    
    
    return_list[[1]] <- geneset_onehot %>% as.list()
    return_list[[2]]<- colnames(geneset_onehot)
    
    return(return_list)
  } ## 1 = geneset_list, 2 = col_name
  fix_load_SKAT <- function(data_name){
    print(paste0(data_name, " fix load!"))
    
    if(data_name == "IPDGC"){
      test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
        as_tibble() %>% select(-CADD13_RawScore, -CADD13_PHRED)
      
      # include indel
      path <- list.files(path = "/home/jinoo/skat-o/SKAT_data/indel_cadd/", full.names = T)
      indel <- lapply(path, FUN = function(x){
        temp <- fread(file = x, header = T, sep = "\t", stringsAsFactors = F) %>% 
          mutate_at(1, funs(as.character(.)))
      }) %>% bind_rows()
      
      indel <- indel %>% rename(CHROM = `#CHROM`, CADD13_RawScore = RawScore,CADD13_PHRED = PHRED)
      
      test_fix <- left_join(x = test_fix, y = indel, by = c("CHROM","POS","REF","ALT")) 
      
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
      test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
      
      test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
      for(change in 1:length(test_id)){
        if(test_id[change] == ""){
          test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
        }
      }
      test_fix$ID <- test_id
      
      freq <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,"_freq.frq"), header = T) %>%
        rename(ID = SNP)
      test_fix <- left_join(x = test_fix, y = freq, by = "ID")
      
      
      
    } else{
      test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
        as_tibble()
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
      test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
      
      test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
      for(change in 1:length(test_id)){
        if(test_id[change] == ""){
          test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
        }
      }
      test_fix$ID <- test_id
      
      freq <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,"_freq.frq"), header = T) %>%
        rename(ID = SNP)
      test_fix <- left_join(x = test_fix, y = freq, by = "ID")
      
    }
    
    # O2 HGNC symobl 0719
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
    
    # CHD HGNC
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^LPHN3$", replacement = "ADGRL3")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^ALS2CR11$", replacement = "C2CD6")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1737$", replacement = "CIPC")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^CCDC129$", replacement = "ITPRID1")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^SUV420H1$", replacement = "KMT5B")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^WIBG$", replacement = "PYM1")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^FAM65C$", replacement = "RIPOR3")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1468$", replacement = "RELCH")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA2018$", replacement = "USF3")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA0196$", replacement = "WASHC5")
    
    {  # ADGRL3  -> LPHN3
      # ALS2CR11 -> C2CD6
      # CIPC -> KIAA1737
      # CNTF -> not fix
      # ITPRID1 -> CCDC129
      # KMT5B -> SUV420H1
      # NUP63 -> not fix
      # PYM1 -> WIBG
      # RELCH -> KIAA1468
      # RIPOR3 -> FAM65C
      # USF3 -> KIAA2018
      # WASHC5 -> KIAA0196
    }
    
    
    return(test_fix)
  } # data_name, "IPDGC", "NeuroX"
  geneset_setid <- function(geneset_merge, col_name, test_fix, data_name, index_, CADD_score){
    print("geneset SetID making!!")
    pb <- progress_bar$new(total = length(index_), clear = F)
    pb$tick(0)
    Sys.sleep(0.01)
    
    for(gene in index_){ # 2:length(col_name), c(2,3,7,16,17)
      geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
      variant <- list()
      ## Parkinson geneset variant
      for(i in 1:length(geneset)){
        variant[[i]] <- subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")),
                               select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED"))
      }
      variant <- bind_rows(variant)
      
      ### 1. nonsynonymous geneset
      nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                 | ExonicFunc.knownGene == "stopgain"
                                                 | ExonicFunc.knownGene ==  "stoploss"
                                                 | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                 | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                 | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                 | ExonicFunc.knownGene ==  "splicing"
      ))
      nonsynonymous <- subset(nonsynonymous, select = "ID")[,1]
      
      
      ### 2. CADD > 12.37 variant
      cadd <- subset(variant, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                         | ExonicFunc.knownGene == "stopgain"
                                         | ExonicFunc.knownGene ==  "stoploss"
                                         | ExonicFunc.knownGene ==  "frameshift_deletion"
                                         | ExonicFunc.knownGene ==  "frameshift_insertion"
                                         | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                         | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED >= CADD_score)
      cadd <- subset(cadd, select = "ID")[,1]
      
      ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
      lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                        | ExonicFunc.knownGene ==  "stoploss"
                                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                        | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED >= CADD_score)
      
      lof <- subset(lof, select = "ID")[,1]
      
      type <- list(nonsynonymous = nonsynonymous, cadd = cadd, lof = lof)
      
      setID <- list()
      for(j in 1:length(type)){
        if(nrow(type[[j]]) == 0){
          next
        }
        temp <- data.frame(TYPE=rep(names(type[j]), length(type[[j]])), stringsAsFactors = F)
        setID[[j]] <- cbind(temp, ID=type[[j]])
      }
      setID <- bind_rows(setID)
      setID$TYPE <- paste0(setID$TYPE, "__",col_name[gene])
      # system("rm -rf skat.SetID")
      fwrite(x = setID, file = paste0(data_name, "_skat.SetID"), sep = "\t",row.names = F, quote = F, col.names = F,append = T)
      pb$tick(1)
      Sys.sleep(0.01)
    }
    
  } 
  gene_setid <- function(index = "default_all", geneset, test_fix, data_name, CADD_score){
    print("gene SetID making!!")
    if(is.character(index)){ # if 1, all gene
      variant_gene <- unique(geneset[[1]]$GENE[!is.na(geneset[[1]]$GENE)])
      
      pb <- progress_bar$new(total = length(variant_gene), clear = F)
      pb$tick(0)
      Sys.sleep(0.01)
      
      #LOF adding
      temp <- list()
      for(i in 1:length(variant_gene)){
        t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED >= CADD_score))
        temp[[i]] <- subset.data.frame(t1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                     | ExonicFunc.knownGene == "stopgain"
                                                     | ExonicFunc.knownGene ==  "stoploss"
                                                     | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                     | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                     | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                     | ExonicFunc.knownGene ==  "splicing"), 
                                       select = c("Gene.knownGene","ID","ID2", "CADD13_PHRED"))
        pb$tick(1)
        Sys.sleep(0.01)
      }
      setID <- bind_rows(temp)
      
      write.table(x = select(setID, Gene.knownGene, ID), file = paste0(data_name, "_skat_gene.SetID"), sep = "\t",row.names = F, quote = F, col.names = F,append = T)
      # write.table(x = select(setID, Gene.knownGene, ID2, CADD13_PHRED), file = paste0(data_name, "_skat_gene_variant.txt"), 
      #             sep = "\t",row.names = F, quote = F, col.names = F,append = T)
    } else{
      
      target_gene <- c()
      
      for(gene in index){
        target_gene <- c(target_gene, geneset[[1]][[gene]][!is.na(geneset[[1]][[gene]])])
      }
      
      variant_gene <- unique.default(target_gene)
      
      pb <- progress_bar$new(total = length(variant_gene), clear = F)
      pb$tick(0)
      Sys.sleep(0.01)
      
      #LOF adding
      temp <- list()
      for(i in 1:length(variant_gene)){
        t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED >= CADD_score))
        temp[[i]] <- subset.data.frame(t1, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                      | ExonicFunc.knownGene == "stopgain"
                                                      | ExonicFunc.knownGene ==  "stoploss"
                                                      | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                      | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                      | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                      | ExonicFunc.knownGene ==  "splicing"), 
                                       select = c("Gene.knownGene","ID","ID2", "CADD13_PHRED"))
        pb$tick(1)
        Sys.sleep(0.01)
      }
      setID <- bind_rows(temp)
      
      write.table(x = select(setID, Gene.knownGene, ID), file = paste0(data_name, "_skat_gene.SetID"), sep = "\t",row.names = F, quote = F, col.names = F,append = T)
      # write.table(x = select(setID, Gene.knownGene, ID2, CADD13_PHRED), file = paste0(data_name, "_skat_gene_variant.txt"), 
      #             sep = "\t",row.names = F, quote = F, col.names = F,append = T)
      
    }
    
  } # all_gene, column 1 of geneset
  run_skat_all_cov <- function(data_name, flag = "geneset", re = 0, add_name){
    print("SKAT run")
    if(flag == "geneset"){
      
      Generate_SSD_SetID(File.Bed = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bed"),
                         File.Bim = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bim"), 
                         File.Fam = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                         File.SetID = paste0(data_name,"_skat.SetID"), 
                         File.SSD = paste0(data_name,".SSD"), 
                         File.Info = paste0(data_name,".INFO"))
      
      
      FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                                File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = FALSE )
      SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
      if(data_name == "IPDGC"){
        # obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + 
        #                         F_MISS, data = FAM, out_type="C", Adjustment = F, n.Resampling = re)
        
        # 20 MDS
        obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                                C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                                C19 + C20 + F_MISS, data = FAM, out_type="C", Adjustment = F, n.Resampling = re)
        
      }else{
        obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + F_MISS,
                              data = FAM, out_type="C", Adjustment = F, n.Resampling = re)
      }
      # skat run
      {
        cl <- makeCluster(detectCores() - 1)
        clusterExport(cl, c("obj","data_name","re"), envir=environment())
        clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
      }
      
      SKAT_result <- 
        parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index){
          FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                                    File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = FALSE)
          SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
          
          Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
          out_05 <- SKAT(Z, obj = obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.05)
          out_03 <- SKAT(Z, obj = obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03)
          out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01)
          
          if(re == 0){
            skat_005 <- tibble(SetID_005 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_05$param$n.marker,
                               n_marker_test = out_05$param$n.marker.test, p_value = out_05$p.value)
            skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker,
                               n_marker_test = out_03$param$n.marker.test, p_value = out_03$p.value)
            skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker,
                               n_marker_test = out_01$param$n.marker.test, p_value = out_01$p.value)
          }else{
            resample_p_05 <- Get_Resampling_Pvalue(out_05)[[1]]
            resample_p_03 <- Get_Resampling_Pvalue(out_03)[[1]]
            resample_p_01 <- Get_Resampling_Pvalue(out_01)[[1]]
            
            skat_005 <- tibble(SetID_005 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_05$param$n.marker, 
                               n_marker_test = out_05$param$n.marker.test, 
                               emprical_pvalue = resample_p_05, p_value = out_05$p.value)
            skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker, 
                               n_marker_test = out_03$param$n.marker.test, 
                               emprical_pvalue = resample_p_03, p_value = out_03$p.value)
            skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker, 
                               n_marker_test = out_01$param$n.marker.test, 
                               emprical_pvalue = resample_p_01, p_value = out_01$p.value)
          }
          
          result <- bind_cols(skat_001, skat_003, skat_005)
          return(result)
          
        })
      
      stopCluster(cl)
      
      results <- SKAT_result %>% bind_rows()
      
      fwrite(x = results, file = paste0(data_name,"_result_cov_all_", add_name,".txt"), col.names = T, row.names = F, sep = "\t")
      file.remove(list.files()[!str_detect(list.files(), pattern = "cov|MAC|SetID")])
      Close_SSD()
      
    } else {
      Generate_SSD_SetID(File.Bed = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bed"),
                         File.Bim = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bim"), 
                         File.Fam = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                         File.SetID = paste0(data_name,"_skat_gene.SetID"), 
                         File.SSD = paste0(data_name,"_gene.SSD"), 
                         File.Info = paste0(data_name,"_gene.INFO"))
      
      
      FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                                File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = FALSE )
      SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
      
      if(data_name == "IPDGC"){
        obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                                C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                                C19 + C20 + F_MISS, data = FAM, out_type="C", Adjustment = T, n.Resampling = re)
      }else{
        obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + F_MISS,
                              data = FAM, out_type="C", Adjustment = F, n.Resampling = re)
      }
      # skat run
      {
        cl <- makeCluster(detectCores() - 1)
        clusterExport(cl, c("re","data_name", "obj"), envir=environment())
        clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
      }
      
      SKAT_result <- 
        parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index){
          FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                                    File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = F)
          SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
          
          Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T);Z[Z==9] <- NA
          
          out_05 <- SKAT(Z, obj = obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.05, estimate_MAF = 2) 
          out_03 <- SKAT(Z, obj = obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03, estimate_MAF = 2)
          out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01, estimate_MAF = 2)
          
          if(re == 0){
            skat_005 <- tibble(SetID_005 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_05$param$n.marker,
                               n_marker_test = out_05$param$n.marker.test, p_value = out_05$p.value)
            skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker,
                               n_marker_test = out_03$param$n.marker.test, p_value = out_03$p.value)
            skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker,
                               n_marker_test = out_01$param$n.marker.test, p_value = out_01$p.value)
          }else{
            resample_p_05 <- Get_Resampling_Pvalue(out_05)[[1]]
            resample_p_03 <- Get_Resampling_Pvalue(out_03)[[1]]
            resample_p_01 <- Get_Resampling_Pvalue(out_01)[[1]]
            
            skat_005 <- tibble(SetID_005 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_05$param$n.marker, 
                               n_marker_test = out_05$param$n.marker.test, 
                               emprical_pvalue = resample_p_05, p_value = out_05$p.value)
            skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker, 
                               n_marker_test = out_03$param$n.marker.test, 
                               emprical_pvalue = resample_p_03, p_value = out_03$p.value)
            skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker, 
                               n_marker_test = out_01$param$n.marker.test, 
                               emprical_pvalue = resample_p_01, p_value = out_01$p.value)
          }
          
          
          result <- bind_cols(skat_001, skat_003, skat_005)
          return(bind_rows(result))
          
        })
      
      
      stopCluster(cl)
      
      SKAT_result %>% bind_rows() %>% 
        fwrite(file = paste0(data_name,"_result_cov_all_gene_",add_name,".txt"), col.names = T, row.names = F, sep = "\t")
      file.remove(list.files()[!str_detect(list.files(), pattern = "cov|MAC|SetID")])
      Close_SSD()
    }
    
  }
  
  MAC_calculation <- function(data_name, add_name){
    print("MAC calculation run")
    Generate_SSD_SetID(File.Bed = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bed"),
                       File.Bim = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".bim"), 
                       File.Fam = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                       File.SetID = paste0(data_name,"_skat_gene.SetID"), 
                       File.SSD = paste0(data_name,"_gene.SSD"), 
                       File.Info = paste0(data_name,"_gene.INFO"))
    
    FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                              File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
    
    {
      cl <- makeCluster(detectCores() - 1)
      clusterExport(cl, c("data_name"), envir=environment())
      clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
    }
    
    MAC_result <- parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index){
      
      FAM <- Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                                File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = F)
      SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
      
      Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T);Z[Z==9] <- NA
      
      #### MAC report
      Z_tibble <- as_tibble(Z) %>% select_if(colMeans(., na.rm = T)/2 != 0 & colMeans(., na.rm = T)/2 <= 0.03)
      control <- which(FAM$Phenotype == 1);case <- which(FAM$Phenotype == 2)
      Z_control <- Z_tibble %>% slice(control);Z_case <- Z_tibble %>% slice(case)
      control_MAC <- Z_control %>% apply(X = ., MARGIN = 1, function(temp) sum(temp, na.rm = T)) %>% sum(na.rm = T)
      case_MAC <- Z_case %>% apply(X = ., MARGIN = 1, function(temp) sum(temp, na.rm = T)) %>% sum(na.rm = T)
      IID_number <- Z_tibble %>% apply(X = ., MARGIN = 1, function(temp) sum(temp, na.rm = T)) %>% 
        .[.>=1] %>% length()
      
      
      result_MAC <- tibble(Gene = SSD.INFO$SetInfo$SetID[SET_index], Variants = ncol(Z_tibble), total_MAC = (control_MAC + case_MAC),
                           N = IID_number, Case_MAC = case_MAC, Control_MAC = control_MAC)
      
      return(result_MAC)
      
    })
    
    stopCluster(cl)
    
    MAC_result %>% bind_rows() %>% 
      write_delim(path = paste0(data_name, "_MAC_gene_", add_name,".txt"), delim = "\t")
    file.remove(list.files()[!str_detect(list.files(), pattern = "cov|MAC|SetID")])
    Close_SSD()
  }
}

## SNP TABLE FUCNTION
{
  # preprocessing
  geneset_load <- function(){
    print("Geneset load_0428!")
    return_list <- list()
    geneset_merge <- fread(file = "/home/jinoo/skat-o/SKAT_data/1 gene sets 190428_exclude_ARX_for SKAT_O.txt", header = T, 
                           sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
    
    one_hot <- lapply(X=select(geneset_merge, -HGNC, -M2, -'M2-GBA-LRRK2', -'M2_Brain_gene', -'M2_NOT_brain'), 
                      FUN = function(x){ifelse(nchar(x) > 0, 1, 0)}) %>% bind_cols()
    geneset_onehot <- cbind(GENE=geneset_merge$HGNC, one_hot) %>% mutate(GENE = as.character(GENE)) %>% 
      cbind(M2 = as.numeric(geneset_merge$M2)) %>% cbind('M2-GBA-LRRK2'= as.numeric(geneset_merge$`M2-GBA-LRRK2`)) %>%
      cbind('M2_Brain_gene' = as.numeric(geneset_merge$`M2_Brain_gene`)) %>% 
      cbind('M2_NOT_brain'= as.numeric(geneset_merge$`M2_NOT_brain`))
    
    geneset_onehot <- geneset_onehot %>% mutate(., 'M2_Brain_gene-LRRK2-GBA' = M2_Brain_gene)
    geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[554]] <- 0
    geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[755]] <- 0
    
    for(geneset in 2:ncol(geneset_onehot)){
      for(k in 1:nrow(geneset_onehot)){
        if(geneset_onehot[k, geneset] >= 1){
          geneset_onehot[k,geneset] <- geneset_onehot[k,1]
        } else{
          geneset_onehot[k,geneset] <- NA
        }
      }
    }
    
    # CHD, DEG adding
    
    CHD <- fread(file = "/home/jinoo/skat-o/SKAT_data/CHD_0424.txt", header = T, 
                 sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble() %>% 
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.)), CHD = NA))
    geneset_onehot <- bind_cols(geneset_onehot, CHD)
    
    DEG <- fread(file = "/home/jinoo/skat-o/SKAT_data/DEG_geneset.txt", header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    meta_DEG <- fread(file = "/home/jinoo/skat-o/SKAT_data/$Table 3 and 4. result_down_intersect_190520_0603.txt", header = T, data.table = F)
    
    
    DEG1 <- DEG$DEG1 %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG1) <- "DEG1"
    DEG2 <- DEG$DEG1YJK[DEG$DEG1YJK != ""] %>% as_tibble() %>%
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG2) <- "DEG1YJK"
    DEG3 <- DEG$DEG1Mito[DEG$DEG1Mito != ""] %>% as_tibble() %>%
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG3) <- "DEG1Mito"
    meta_DEG <- meta_DEG$SYMBOL[meta_DEG$SYMBOL != ""] %>% as_tibble() %>%
      bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(meta_DEG) <- "META_DEG"
    
    geneset_onehot <- geneset_onehot %>% bind_cols(., DEG1) %>% bind_cols(., DEG2) %>% bind_cols(., DEG3) %>%
      bind_cols(., meta_DEG)
    
    
    return_list[[1]] <- geneset_onehot %>% as.list()
    return_list[[2]]<- colnames(geneset_onehot)
    
    return(return_list)
  }
  fix_load <- function(data_name, type){
    print(paste0(data_name, " fix load!"))
    
    
    if(type == "ROW"){
      test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_row_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
        as_tibble() 
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
      
      # O2 HGNC symobl 0719
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
      
      # CHD HGNC
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^LPHN3$", replacement = "ADGRL3")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^ALS2CR11$", replacement = "C2CD6")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1737$", replacement = "CIPC")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^CCDC129$", replacement = "ITPRID1")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^SUV420H1$", replacement = "KMT5B")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^WIBG$", replacement = "PYM1")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1468$", replacement = "RELCH")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA2018$", replacement = "USF3")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA0196$", replacement = "WASHC5")
      test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^FAM65C$", replacement = "RIPOR3")
      
      return(test_fix)
    }
    
    
    if(data_name == "IPDGC" & type != "ROW"){
      test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
        as_tibble() %>% select(-CADD13_RawScore, -CADD13_PHRED)
      
      # include indel
      path <- list.files(path = "/home/jinoo/skat-o/SKAT_data/indel_cadd/", full.names = T)
      indel <- lapply(path, FUN = function(x){
        temp <- fread(file = x, header = T, sep = "\t", stringsAsFactors = F) %>% 
          mutate_at(1, funs(as.character(.)))
      }) %>% bind_rows()
      
      indel <- indel %>% rename(CHROM = `#CHROM`, CADD13_RawScore = RawScore,CADD13_PHRED = PHRED)
      
      test_fix <- left_join(x = test_fix, y = indel, by = c("CHROM","POS","REF","ALT")) 
      
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
      test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
      
      test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
      for(change in 1:length(test_id)){
        if(test_id[change] == ""){
          test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
        }
      }
      test_fix$ID <- test_id
      
      freq <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,"_freq.frq"), header = T) %>%
        rename(ID = SNP)
      test_fix <- left_join(x = test_fix, y = freq, by = "ID")
      
    } else if(data_name == "NeuroX" & type != "ROW"){
      test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
        as_tibble()
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
      test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
      test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
      
      test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
      for(change in 1:length(test_id)){
        if(test_id[change] == ""){
          test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
        }
      }
      test_fix$ID <- test_id
      
      freq <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,"_freq.frq"), header = T) %>%
        rename(ID = SNP)
      test_fix <- left_join(x = test_fix, y = freq, by = "ID")
      
      
    }
    
    # O2 HGNC symobl 0719
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
    
    # CHD HGNC
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^LPHN3$", replacement = "ADGRL3")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^ALS2CR11$", replacement = "C2CD6")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1737$", replacement = "CIPC")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^CCDC129$", replacement = "ITPRID1")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^SUV420H1$", replacement = "KMT5B")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^WIBG$", replacement = "PYM1")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA1468$", replacement = "RELCH")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA2018$", replacement = "USF3")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^KIAA0196$", replacement = "WASHC5")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^FAM65C$", replacement = "RIPOR3")
    
    return(test_fix)
  } # data_name, "IPDGC", "NeuroX"
  geneset_extract <- function(geneset_merge, col_name, test_fix, data_name, index_){
    print("geneset SetID making!!")
    
    for(gene in index_){ # 2:length(col_name), c(2,3,7,16,17)
      geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
      variant <- list()
      ## Parkinson geneset variant
      for(i in 1:length(geneset)){
        variant[[i]] <- subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")),
                               select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED"))
      }
      variant <- bind_rows(variant)
      
      ### 1. nonsynonymous geneset
      nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"))
      fwrite(x = subset(nonsynonymous, select = "ID")[,1], file = paste0(data_name,"_",col_name[gene],"_non-syn.txt"), col.names = F)
      
      
      
      ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
      lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                        | ExonicFunc.knownGene ==  "stoploss"
                                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                        | ExonicFunc.knownGene ==  "splicing"))
      
      fwrite(x = subset(lof, select = "ID")[,1], file = paste0(data_name,"_",col_name[gene],"_lof.txt"), col.names = F)
    }
    
  } 
  clinvar_load <- function(){
    # clinvar preprocessing
    clinvar <- fread(file = "/home/jinoo/skat-o/SKAT_data/clinvar_20190506.vcf", header = T) %>%
      mutate(., CHROM = str_replace(string = CHROM, pattern = "X", replacement = "23")) %>%
      mutate(., CHROM = str_replace(string = CHROM, pattern = "Y", replacement = "24")) %>%
      mutate(., CHROM = str_replace(string = CHROM, pattern = "MT", replacement = "26")) %>%
      mutate(., CHROMPOS = paste(CHROM, POS, sep = ":")) %>%
      mutate(., CHROMPOS_R_A = paste0(CHROMPOS, REF, ALT))
    
    info <- mclapply(X = clinvar$INFO, FUN = clinvar_info, mc.cores = detectCores() - 1) %>% bind_rows()
    clinvar <- clinvar %>% select(-INFO) %>% bind_cols(info)
    return(clinvar)
  }
  clinvar_info <- function(temp){
    str_split(temp,pattern = ";")[[1]] %>% 
      as_tibble() %>% separate(value, c("info","value"), sep = "=") %>% 
      spread(key = "info", value = "value") %>% return()
  }
  re_CLNSIG <- function(CHROMPOS, clinvar){
    select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], CLNSIG)[,1] %>%
      return()
  }
  re_CLNSIGCONF <- function(CHROMPOS, clinvar){
    select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], CLNSIGCONF)[,1] %>%
      return()
  }
  
  table3_CLSIG <- function(table){
    result <- mclapply(X = unique.default(table$IID), function(target_IID){
      select_IID <- table %>% filter(IID == target_IID) %>% select(IID, PHENOTYPE) %>% distinct()
      
      clsig_homo <- NULL
      # zigosity, homo count * 2, hetero * 1
      if(nrow(filter(table, IID == target_IID, Homozygous_M == T)) > 0){
        clsig_homo <- table %>% filter(IID == target_IID, Homozygous_M == T) %>%  
          count(Clinical_Significance) %>% spread(key = "Clinical_Significance", value = "n") %>%
          apply(., MARGIN = 2, FUN = function(x) x * 2)
      }
      clsig_hetero <- table %>% filter(IID == target_IID, Heterozygous_M == T) %>% 
        count(Clinical_Significance) %>% spread(key = "Clinical_Significance", value = "n")  
      
      clsig <- bind_rows(clsig_homo, clsig_hetero);clsig[is.na(clsig)] <- 0
      clsig <- clsig %>% apply(., MARGIN = 2, sum) %>% t() %>% as_tibble()
      
      return(bind_cols(select_IID, clsig))
    }, mc.cores = detectCores() - 1)
    return(result)  
  }
  table3_calc <- function(table_list){
    result_list <- list()
    col_name <- colnames(table_list)
    
    # if(nrow(table_list) == 0) return(tibble())
    
    for(index in 3:ncol(table_list)){
      temp <- table_list %>% select(col_name[index]) %>% 
        group_by_all() %>% 
        summarise(value = n()) %>%
        mutate_at(1, funs(as.character(.)))
      
      if(nrow(temp) == 0){
        result <- tibble(rep(0,9), .name_repair = ~col_name[index])
      } else{
        zero_row <- temp[,1] %in% "0"
        if(zero_row[1] == F)
          temp <- bind_rows(tibble("0", 0, .name_repair = ~c(col_name[index], "value")), temp)
        
        if(nrow(temp) <= 8){
          temp <- tibble(c(as.character(0:7), "8+"), .name_repair = ~c(col_name[index])) %>% 
            left_join(x = ., y = temp, by = col_name[index])
        } else if(nrow(temp) >= 9){
          more_value <- temp %>% slice(9:nrow(.)) %>%
            pull(2) %>% sum()
          temp <- bind_rows(slice(temp, 1:8), tibble("8+", more_value, .name_repair = ~ c(col_name[index], "value")))
        }
        
        temp$value[is.na(temp$value)] <- 0
        result <- temp[2];colnames(result) <- col_name[index]
      }
      
      result_list[[index-2]] <- result
    }
    return(result_list)
  }
  
  # result processing
  snp_table <- function(data_name = "IPDGC", index_, clinvar){
    core <- detectCores() - 1
    geneset_result <- list()
    geneset <- geneset_load()
    fix <- fix_load(data_name, type = "QC")
    
    geneset_extract(geneset[[1]], geneset[[2]], fix, data_name, index_)
    
    for(index in index_){
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_non-syn.txt --recode A --out test_dosage_nonsyn", 
                  data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_lof.txt --recode A --out test_dosage_lof", 
                  data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
      
      nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, SEX) %>% as_tibble()
      lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, SEX) %>% as_tibble()
      
      nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
      # nonsyn_dosage[which(nonsyn_dosage$IID == "UMARY-1571"), which(colnames(nonsyn_dosage) == "rs138008832_A")] <- 1 ######### remove
      lof_dosage[is.na(lof_dosage)] <- 0
      
      print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
      nonsyn_result <- mclapply(X = 4:ncol(nonsyn_dosage), FUN = function(col_len){
        anno_data_nonsyn <- filter(fix, ID == str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]) %>%
          select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
        
        sample_nonsyn <- select(nonsyn_dosage, c(1:3, col_len)) %>% 
          # filter(., .[,3] >= 1) %>% 
          mutate(Heterozygous_M = ifelse(.[,4] == 1, TRUE, FALSE),
                 Homozygous_M = ifelse(.[,4] == 2, TRUE, FALSE),
                 Heterozygous = ifelse(.[,4] == 0, TRUE, FALSE)) %>%
          select(., IID, PHENOTYPE, SEX, Heterozygous , Heterozygous_M, Homozygous_M)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_nonsyn) == 0) return(anno_data_mul)
        for(i in 1:nrow(sample_nonsyn))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_nonsyn)  
        return(bind_cols(sample_nonsyn, anno_data_mul))
      }, mc.cores = core) %>% bind_rows()
      
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 4:ncol(lof_dosage), FUN = function(col_len){
        anno_data_lof <- filter(fix, ID == str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]) %>%
          select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
        
        sample_lof <- select(lof_dosage, c(1:3, col_len)) %>% 
          # filter(., .[,3] >= 1) %>% 
          mutate(Heterozygous_M = ifelse(.[,4] == 1, TRUE, FALSE),
                 Homozygous_M = ifelse(.[,4] == 2, TRUE, FALSE),
                 Heterozygous = ifelse(.[,4] == 0, TRUE, FALSE)) %>%
          select(., IID, PHENOTYPE, SEX, Heterozygous, Heterozygous_M, Homozygous_M)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_lof) == 0) return(anno_data_mul)
        
        for(i in 1:nrow(sample_lof))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_lof)    
        return(bind_cols(sample_lof, anno_data_mul))
      }, mc.cores = core) %>% bind_rows()
      
      nonsyn_result <- nonsyn_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               Datasets = data_name, geneset = geneset[[2]][index], 
               Func = "Nonsynonymous")
      
      lof_result <- lof_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               Datasets = data_name, 
               geneset = geneset[[2]][index], Func = "LoF")
      
      result <- bind_rows(nonsyn_result, lof_result)
      
      result_DF <- mclapply(X = 1:nrow(result), FUN = function(rows){
        temp <- result[rows, ]
        if(temp$SEX == 1 & (temp$CHROM == "23" | temp$CHROM == "24") & temp$Homozygous_M == TRUE){
          temp$Heterozygous_M <- TRUE
          temp$Homozygous_M <- FALSE
        }
        
        return(temp)
      }, mc.cores = detectCores() -1) %>% bind_rows()
      
      geneset_result[[geneset[[2]][index]]] <- result_DF
    }
    
    system("rm -rf *.log");system("rm -rf *.raw");system("rm -rf *.hh")
    system("rm -rf *.txt")
    return(geneset_result)
  } # gender adding 0722
  snp_table_ALL <- function(data_name = "IPDGC", index_, clinvar){
    core <- detectCores() - 1
    geneset_result <- list()
    geneset <- geneset_load()
    fix <- fix_load(data_name, type = "QC")
    
    geneset_extract(geneset[[1]], geneset[[2]], fix, data_name, index_)
    
    if(data_name == "IPDGC"){
      for(index in index_){
        system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_non-syn.txt --recode A --out test_dosage_nonsyn", 
                    data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
        system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_lof.txt --recode A --out test_dosage_lof", 
                    data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
        
        nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
        lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
        
        nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
        lof_dosage[is.na(lof_dosage)] <- 0
        
        print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
        nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
          anno_data_nonsyn <- filter(fix, ID == str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]) %>%
            select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
            mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                   CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>%
          bind_rows()
        
        
        print(paste(geneset[[2]][index],"lof", sep = " "))
        lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
          anno_data_lof <- filter(fix, ID == str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]) %>%
            select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
            mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                   CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>% 
          bind_rows()
        
        nonsyn_result <- nonsyn_result %>%
          mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
                 geneset = geneset[[2]][index], Func = "Nonsynonymous")
        
        lof_result <- lof_result %>%
          mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT)%in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
                 geneset = geneset[[2]][index], Func = "LoF")
        geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
      }
      
    } else{
      
      for(index in index_){
        system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_non-syn.txt --recode A --out test_dosage_nonsyn", 
                    data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
        system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_lof.txt --recode A --out test_dosage_lof", 
                    data_name = data_name, geneset = geneset[[2]][index]), ignore.stdout = T)
        
        nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
        lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
        
        nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
        lof_dosage[is.na(lof_dosage)] <- 0
        
        print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
        nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
          anno_data_nonsyn <- filter(fix, ID == str_sub(colnames(nonsyn_dosage)[col_len], end = -3)) %>%
            select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
            mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                   CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>% 
          bind_rows()
        
        
        print(paste(geneset[[2]][index],"lof", sep = " "))
        lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
          anno_data_lof <- filter(fix, ID == str_sub(colnames(lof_dosage)[col_len], end = -3)) %>%
            select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
            mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                   CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>% 
          bind_rows()
        
        nonsyn_result <- nonsyn_result %>%
          mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
                 geneset = geneset[[2]][index], Func = "Nonsynonymous")
        
        lof_result <- lof_result %>%
          mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
                 geneset = geneset[[2]][index], Func = "LoF")
        geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
      }
      
    }
    
    system("rm -rf *.log");system("rm -rf *.raw");system("rm -rf *.hh")
    system("rm -rf *.txt")
    return(geneset_result)
  }
  
  table3_ver_1 <- function(data_list, type = "Pathogenic", geneset_name, MAF_value = 0.03, CADD_score = 20){
    index <- which(names(data_list) == geneset_name)
    
    temp <- data_list[[index]]
    temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                              pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                              replacement = "Pathogenic_ALL")
    temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
    
    if(type == "Pathogenic"){
      temp2 <- temp %>% arrange(IID) %>% 
        filter(MAF <= MAF_value, Heterozygous != TRUE, Clinical_Significance == "Pathogenic_ALL") %>%
        table3_CLSIG() %>% bind_rows()
      temp2[is.na(temp2)] <- 0
      
      control <- temp2 %>% filter(PHENOTYPE == 1) %>%
        rename(Pathogenic_ALL_control = Pathogenic_ALL) %>%
        table3_calc() 
      case <- temp2 %>% filter(PHENOTYPE == 2) %>%
        rename(Pathogenic_ALL_case = Pathogenic_ALL) %>%
        table3_calc() 
      
    } else if(type == "LoF"){
      temp2 <- temp %>% arrange(IID) %>% 
        filter(MAF <= MAF_value, Func == "LoF", Heterozygous != TRUE, Clinical_Significance != "Pathogenic_ALL") %>%
        table3_CLSIG() %>% bind_rows()
      temp2[is.na(temp2)] <- 0
      
      lof_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(lof_sum) <- "LoF"
      temp2 <- bind_cols(temp2[, 1:2], lof_sum)
      
      control <- temp2 %>% filter(PHENOTYPE == 1) %>%
        rename(LoF_control = LoF) %>%
        table3_calc() 
      case <- temp2 %>% filter(PHENOTYPE == 2) %>%
        rename(LoF_case = LoF) %>%
        table3_calc() 
      
    }else if(type == "CADD_MAF"){
      temp2 <- temp %>% arrange(IID) %>% 
        filter(MAF <= MAF_value, Heterozygous != TRUE, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL") %>%
        table3_CLSIG() %>% bind_rows()
      temp2[is.na(temp2)] <- 0
      
      cadd_maf_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(cadd_maf_sum) <- "cadd_maf"
      temp2 <- bind_cols(temp2[, 1:2], cadd_maf_sum)
      
      control <- temp2 %>% filter(PHENOTYPE == 1) %>%
        rename(CADD_MAF_control = cadd_maf) %>%
        table3_calc() 
      
      case <- temp2 %>% filter(PHENOTYPE == 2) %>%
        rename(CADD_MAF_case = cadd_maf) %>%
        table3_calc() 
      
      
      
    } 
    
    
    result <- bind_cols(case, control) 
    
    if(result[1,1] == 0)
      result[1,1] <- 443 - result[,1] %>% sum
    if(result[1,2] == 0)
      result[1,2] <- 333 - result[,2] %>% sum
    
    return(result)
  }
  table3_ver_2_gene <- function(WES_table, geneset_name, CADD_score = 20){
    index <- which(names(WES_table) == geneset_name)
    
    temp <- WES_table[[index]]
    temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                              pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                              replacement = "Pathogenic_ALL")
    temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
    gene <- temp$Gene.knownGene %>% unique() %>% sort()
    result <- list()
    category <- geneset_name
    # Genes in ~ 
    
    for(Genes in gene){
      Pathogenic_ALL <- temp %>% 
        filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", Gene.knownGene == Genes) %>% select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
      
      LoF_003 <- temp %>% 
        filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", Gene.knownGene == Genes, Func == "LoF") %>% select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      LoF_001 <- temp %>% 
        filter(MAF <= 0.01, Clinical_Significance != "Pathogenic_ALL", Gene.knownGene == Genes, Func == "LoF") %>% select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      
      CADD13_MAF003 <- temp %>% 
        filter(MAF <= 0.03, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", Gene.knownGene == Genes) %>% 
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      CADD13_MAF001 <- temp %>% 
        filter(MAF <= 0.01, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", Gene.knownGene == Genes) %>% 
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      
      result[[`Genes`]] <- tibble(category, Genes, Pathogenic_ALL, LoF_003, CADD13_MAF003, LoF_001, CADD13_MAF001) 
    }
    
    return(result)
  }
  table3_ver_3_sample <- function(data_list, geneset_name, CADD_score = 20){
    index <- which(names(data_list) == geneset_name)
    temp <- data_list[[index]]
    temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                              pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                              replacement = "Pathogenic_ALL")
    temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
    sample <- temp$IID %>% unique() %>% sort()
    
    sample_result <- mclapply(X = sample, FUN = function(iid){
      # result <- list()
      PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
        unique() 
      PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
      SEX <- temp %>% filter(IID == iid) %>% select(SEX) %>% pull(1) %>%
        unique()
      
      # Pathogenic
      {
        Pathogenic_Hetro <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Heterozygous_M == T) %>% 
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
        
        if(Pathogenic_Hetro != 0){
          Hetero_name <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                   IID == iid, Heterozygous_M == T) %>% .$ID
          Pathogenic_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
        } else {Pathogenic_Hetro_ID <- " "}
        
        
        Pathogenic_Homo <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Homozygous_M == T) %>% 
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
        
        if(Pathogenic_Homo != 0){
          Homo_name <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                   IID == iid, Homozygous_M == T) %>% .$ID
          Pathogenic_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
        } else {Pathogenic_Homo_ID <- " "}
        }
      
      # LoF
      {
        LoF_Hetero <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Func == "LoF", Heterozygous_M == T) %>%
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
        
        if(LoF_Hetero != 0){
          Hetero_name <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$ID
          LoF_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
        } else {LoF_Hetro_ID <- " "}
        
        LoF_Homo <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Func == "LoF", Homozygous_M == T) %>%
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
        
        if(LoF_Homo != 0){
          Homo_name <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Func == "LoF", Homozygous_M == T) %>% .$ID
          LoF_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
        } else {LoF_Homo_ID <- " "}
      }
      
      # CADD
      {
        for(value in c(0.03, 0.01)){
          CADD_Hetero <- temp %>% 
            filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Heterozygous_M == T) %>%
            select(Clinical_Significance) %>%
            group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
          
          if(CADD_Hetero != 0){
            Hetero_name <- temp %>% 
              filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                     IID == iid, Heterozygous_M == T) %>% .$ID
            CADD_Hetro_name <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
          } else {CADD_Hetro_name <- " "}
          
          CADD_Homo <- temp %>% 
            filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Homozygous_M == T) %>%
            select(Clinical_Significance) %>%
            group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
          
          if(CADD_Homo != 0){
            Homo_name <- temp %>% 
              filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                     IID == iid, Homozygous_M == T) %>% .$ID
            CADD_Homo_name <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
          } else {CADD_Homo_name <- " "}
          
          if(value == 0.03){
            CADD_MAF003_Hetero <- CADD_Hetero
            CADD_MAF003_Hetero_ID <- CADD_Hetro_name
            CADD_MAF003_Homo <- CADD_Homo
            CADD_MAF003_Homo_ID <- CADD_Homo_name
          } else{
            CADD_MAF001_Hetero <- CADD_Hetero
            CADD_MAF001_Hetero_ID <- CADD_Hetro_name
            CADD_MAF001_Homo <- CADD_Homo
            CADD_MAF001_Homo_ID <- CADD_Homo_name
          }
          
        }
      }
      
      
      tibble(Subject_ID = iid, SEX, PHENOTYPE, Pathogenic_ALL = (Pathogenic_Hetro + (Pathogenic_Homo * 2)),
             LoF_ALL = (LoF_Hetero + (LoF_Homo * 2)),
             CADD_MAF003_ALL = (CADD_MAF003_Hetero + (CADD_MAF003_Homo * 2)),
             CADD_MAF001_ALL = (CADD_MAF001_Hetero + (CADD_MAF001_Homo * 2)),
             Pathogenic_Hetro, Pathogenic_Hetro_ID,
             Pathogenic_Homo, Pathogenic_Homo_ID,
             LoF_Hetero, LoF_Hetro_ID, LoF_Homo, LoF_Homo_ID,
             CADD_MAF003_Hetero, CADD_MAF003_Hetero_ID,
             CADD_MAF003_Homo, CADD_MAF003_Homo_ID,
             CADD_MAF001_Hetero, CADD_MAF001_Hetero_ID,
             CADD_MAF001_Homo, CADD_MAF001_Homo_ID) %>% return()
    }, mc.cores = detectCores() -1)
    
    return(sample_result)
  } 
  table3_ver_3_sample_gene <- function(data_list, geneset_name, type, CADD_score = 20){
    index <- which(names(data_list) == geneset_name)
    temp <- data_list[[index]]
    temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                              pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                              replacement = "Pathogenic_ALL")
    temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
    sample <- temp$IID %>% unique() %>% sort()
    
    sample_result <- mclapply(X = sample, FUN = function(iid){
      # result <- list()
      PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
        unique() 
      PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
      
      # Pathogenic
      
      Pathogenic <- temp %>% 
        filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Heterozygous != T) %>% 
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
      
      if(Pathogenic != 0){
        Pathogenic_name_He <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Heterozygous_M == T) %>% .$ID %>% as_tibble()
        Pathogenic_ID_He <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Heterozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
          bind_cols(., Pathogenic_name_He) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
          t() %>% as_tibble() %>% rename_all(.funs = function(name){
            paste0("Pathogenic_Hetero_",name)
          })
        
        Pathogenic_name_Ho <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Homozygous_M == T) %>% .$ID %>% as_tibble()
        
        if(nrow(Pathogenic_name_Ho) == 0){
          Pathogenic_ID_Ho <- NULL
        } else{
          Pathogenic_ID_Ho <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                   IID == iid, Homozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
            bind_cols(., Pathogenic_name_Ho) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
            t() %>% as_tibble() %>% rename_all(.funs = function(name){
              paste0("Pathogenic_Homo_",name)
            })
        }
        
        Pathogenic_ID <- bind_cols(Pathogenic_ID_He, Pathogenic_ID_Ho)
        
      } else {Pathogenic_ID <- tibble(.rows = 1)}
      
      
      # LoF
      
      LoF <- temp %>% 
        filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
               IID == iid, Func == "LoF", Heterozygous != T) %>%
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      if(LoF!= 0){
        LoF_name_He <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$ID %>% as_tibble()
        
        LoF_ID_He <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
          bind_cols(., LoF_name_He) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
          t() %>% as_tibble() %>% rename_all(.funs = function(name){
            paste0("LoF_Hetero_",name)
          })
        
        LoF_name_Ho <- temp %>% 
          filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Func == "LoF", Homozygous_M == T) %>% .$ID %>% as_tibble()
        
        if(nrow(LoF_name_Ho) == 0){
          LoF_ID_Ho <- NULL
        } else {
          LoF_ID_Ho <- temp %>% 
            filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Func == "LoF", Homozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
            bind_cols(., LoF_name_Ho) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
            t() %>% as_tibble() %>% rename_all(.funs = function(name){
              paste0("LoF_Homo_",name)
            })
        }
        
        LoF_ID <- bind_cols(LoF_ID_He, LoF_ID_Ho)
        
      } else {
        LoF_ID <- tibble(.rows = 1)
      }
      
      
      
      # CADD
      for(value in c(0.03, 0.01)){
        CADD <- temp %>% 
          filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                 IID == iid, Heterozygous != T) %>% select(Clinical_Significance) %>%        
          group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
        
        if(CADD != 0){
          CADD_name_He <- temp %>% 
            filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Heterozygous_M == T) %>% .$ID %>% as_tibble()
          CADD_ID_He <- temp %>% 
            filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Heterozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
            bind_cols(., CADD_name_He) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
            t() %>% as_tibble() %>% rename_all(.funs = function(name){
              paste0("CADD_MAF003_Hetero_",name)
            })
          
          
          CADD_name_Ho <- temp %>% 
            filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                   IID == iid, Homozygous_M == T) %>% .$ID %>% as_tibble()
          
          if(nrow(CADD_name_Ho) == 0){
            CADD_ID_Ho <- NULL
          } else{
            CADD_ID_Ho <- temp %>% 
              filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
                     IID == iid, Homozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
              bind_cols(., CADD_name_Ho) %>% transmute(gene_rs = paste0(value, "(", value1,")")) %>%
              t() %>% as_tibble() %>% rename_all(.funs = function(name){
                paste0("CADD_MAF003_Homo_",name)
              })
          }
          
          CADD_ID <- bind_cols(CADD_ID_He, CADD_ID_Ho)
          
        } else {
          CADD_ID <- tibble(.rows = 1)
        }
        
        
        if(value == 0.03){CADD_MAF003_ID <- CADD_ID
        }else{CADD_MAF001_ID <- CADD_ID}
        
      }
      
      #######
      snp_table <- switch(type, "Pathogenic" = Pathogenic_ID,
                          "LoF" = LoF_ID,
                          "CADD_MAF003" = CADD_MAF003_ID,
                          "CADD_MAF001" = CADD_MAF001_ID,
                          "ALL" = bind_cols(Pathogenic_ID, CADD_MAF003_ID) %>% as_tibble()
      )
      
      tibble(Subject_ID = iid, PHENOTYPE) %>% bind_cols(snp_table) %>% return()
    }, mc.cores = detectCores() -1)
    
    return(sample_result)
  } 
  
  # corr & plot
  {
    plot_data_load <- function(){
      plot_list <- list()
      
      plot_list[["O2"]] <- fread(file = "O2_sample_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>% .[,c(1:6, 23)]
      plot_list[["O2_PARK"]] <- fread(file = "O2_PARK_sample_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>% .[,c(1:6, 23)]
      plot_list[["O2_NONPARK"]] <- fread(file = "O2_NONPARK_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>% .[,c(1:6, 23)]
      
      # O2 %>% filter(PHENOTYPE == "case") %>% base::summary()
      # O2 %>% filter(PHENOTYPE == "control") %>% base::summary()
      
      
      return(plot_list)
    }
    plot_create <- function(data, data_col, geneset){
      
      for(phenotype in c("case","control")){
        plot_df <- filter(data, PHENOTYPE == phenotype)
        
        if(data_col == "Pathogenic_ALL"){
          ggplot(data = plot_df, aes(x = Pathogenic_ALL, y = AGE, group = Pathogenic_ALL)) +
            geom_boxplot() +
            coord_cartesian(xlim = seq(0, 8, 1), ylim = seq(10, 100, 10), expand = T) +
            scale_x_continuous(breaks = seq(0, 8, 1)) +
            scale_y_continuous(breaks = seq(0, 100, 10)) +
            ggtitle(paste0(data_col, "(", phenotype,")","_", geneset)) +
            labs(x = "Variant count", y = "AGE") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20, color = "darkblue", family = "serif"),
                  axis.title = element_text(size = 14, color = "darkblue"),
                  axis.text = element_text(face = "bold", size = 10)) 
          
          ggsave(paste0(data_col, "(", phenotype,")","_", geneset, ".png"))
          
        } else if(data_col == "LoF_ALL"){
          ggplot(data = plot_df, aes(x = LoF_ALL, y = AGE, group = LoF_ALL)) +
            geom_boxplot() +
            coord_cartesian(xlim = seq(0, 8, 1), ylim = seq(10, 100, 10), expand = T) +
            scale_x_continuous(breaks = seq(0, 8, 1)) +
            scale_y_continuous(breaks = seq(0, 100, 10)) +
            ggtitle(paste0(data_col, "(", phenotype,")","_", geneset)) +
            labs(x = "Variant count", y = "AGE") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20, color = "darkblue", family = "serif"),
                  axis.title = element_text(size = 14, color = "darkblue"),
                  axis.text = element_text(face = "bold", size = 10))
          
          ggsave(paste0(data_col, "(", phenotype,")","_", geneset, ".png"))
          
        } else if(data_col == "CADD_MAF003_ALL"){
          ggplot(data = plot_df, aes(x = CADD_MAF003_ALL, y = AGE, group = CADD_MAF003_ALL)) +
            geom_boxplot() +
            coord_cartesian(xlim = seq(0, 8, 1), ylim = seq(10, 100, 10), expand = T) +
            scale_x_continuous(breaks = seq(0, 8, 1)) +
            scale_y_continuous(breaks = seq(0, 100, 10)) +
            ggtitle(paste0(data_col, "(", phenotype,")","_", geneset)) +
            labs(x = "Variant count", y = "AGE") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20, color = "darkblue", family = "serif"),
                  axis.title = element_text(size = 14, color = "darkblue"),
                  axis.text = element_text(face = "bold", size = 10))
          
          ggsave(paste0(data_col, "(", phenotype,")","_", geneset, ".png"))
          
        } else if(data_col == "CADD_MAF001_ALL"){
          
          ggplot(data = plot_df, aes(x = CADD_MAF001_ALL, y = AGE, group = CADD_MAF001_ALL)) +
            geom_boxplot() +
            coord_cartesian(xlim = seq(0, 8, 1), ylim = seq(10, 100, 10), expand = T) +
            scale_x_continuous(breaks = seq(0, 8, 1)) +
            scale_y_continuous(breaks = seq(0, 100, 10)) +
            ggtitle(paste0(data_col, "(", phenotype,")","_", geneset)) +
            labs(x = "Variant count", y = "AGE") +
            theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20, color = "darkblue", family = "serif"),
                  axis.title = element_text(size = 14, color = "darkblue"),
                  axis.text = element_text(face = "bold", size = 10)) 
          ggsave(paste0(data_col, "(", phenotype,")","_", geneset, ".png"))
          
        } else{
          stop("uncorrected name!!!!")
        }
      }
      
    }
    # plot_create(data = O2_NONPARK, data_col = colnames(O2)[3], geneset = "O2_NONPARK")
    plot_result <- function(list_DF, geneset_name){
      
      for(index in 3:(length(list_DF) - 1))
        plot_create(data = list_DF, data_col = colnames(list_DF)[index], geneset = geneset_name)
      
    }
    
    # plot_result(list_DF = plot_list[[1]], geneset_name = "O2")
    # plot_result(list_DF = plot_list[[2]], geneset_name = "O2_PARK")
    # plot_result(list_DF = plot_list[[3]], geneset_name = "O2_NONPARK")
    
    # fread(file = "O2_sample_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>%
    #   select(Subject_ID, PHENOTYPE, AGE,Pathogenic_ALL, LoF_ALL, CADD_MAF003_ALL, CADD_MAF001_ALL) %>% corr_geneset()
    # 
    # fread(file = "O2_PARK_sample_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>%
    #   select(Subject_ID, PHENOTYPE, AGE,Pathogenic_ALL, LoF_ALL, CADD_MAF003_ALL, CADD_MAF001_ALL) %>% corr_geneset()
    # 
    # fread(file = "O2_NONPARK_snv_.txt", header = T, sep = "\t") %>% as_tibble() %>%
    #   select(Subject_ID, PHENOTYPE, AGE,Pathogenic_ALL, LoF_ALL, CADD_MAF003_ALL, CADD_MAF001_ALL) %>% corr_geneset()
    
    corr_geneset <- function(geneset_df){
      return_list <- list()
      return_list[["case"]] <- filter(geneset_df, PHENOTYPE == "case") %>% corr_t()
      return_list[["control"]] <- filter(geneset_df, PHENOTYPE == "control") %>% corr_t()
      
      
      
      return(return_list)
      
    }
    corr_t <- function(table){
      result <- c()
      result <- c(result, cor(x = table$Pathogenic_ALL, y = table$AGE))
      result <- c(result, cor(x = table$LoF_ALL, y = table$AGE))
      result <- c(result, cor(x = table$CADD_MAF003_ALL, y = table$AGE))
      result <- c(result, cor(x = table$CADD_MAF001_ALL, y = table$AGE))
      
      names(result) <- c("Pathogenic", "LoF", "CADD20_MAF003","CADD20_MAF001")
      
      return(result)
    }  
  }
}
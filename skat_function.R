library(SKAT);library(tidyverse);library(data.table)
library(parallel);library(doParallel);library(progress)


# function
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
} ## 1 = geneset_list, 2 = col_name

fix_load <- function(data_name){
  print(paste0(data_name, " fix load!"))
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
  
  return(test_fix)
} # data_name, "IPDGC", "NeuroX"

geneset_setid <- function(geneset_merge, col_name, test_fix, data_name, index_){
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
                                       | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
    cadd <- subset(cadd, select = "ID")[,1]
    
    ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
    lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                      | ExonicFunc.knownGene ==  "stoploss"
                                      | ExonicFunc.knownGene ==  "frameshift_deletion"
                                      | ExonicFunc.knownGene ==  "frameshift_insertion"
                                      | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                      | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
    
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

gene_setid <- function(index = "default_all", geneset, test_fix, data_name){
  print("gene SetID making!!")
  if(is.character(index)){ # if 1, all gene
    variant_gene <- unique(geneset[[1]]$GENE[!is.na(geneset[[1]]$GENE)])
    
    pb <- progress_bar$new(total = length(variant_gene), clear = F)
    pb$tick(0)
    Sys.sleep(0.01)
    
    temp <- list()
    for(i in 1:length(variant_gene)){
      t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED > 12.37))
      temp[[i]] <- subset.data.frame(t1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" ), 
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
    
    temp <- list()
    for(i in 1:length(variant_gene)){
      t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED > 12.37))
      temp[[i]] <- subset.data.frame(t1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" ), 
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

run_skat_all_cov <- function(data_name, flag = "geneset", re = 0){
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
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
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
        FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                                  File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE)
        SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
        
        
        
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
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all_gene.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
  }
  
}

run_skat_all_common_rare_cov <- function(data_name, flag = "geneset", re = 0){
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
        out_common_rare <- SKAT_CommonRare(Z, obj, method = "C", CommonRare_Cutoff = 0.05, r.corr.rare = 1,
                                           r.corr.common = 1, test.type = "Common.Only", missing_cutoff = 0.9)
        
        if(re == 0){
          skat_common_rare <- tibble(SetID_common_rare = SSD.INFO$SetInfo[SET_index,2],
                                     n_marker = out_common_rare$param$n.marker,
                                     n_marker_test = out_common_rare$param$n.marker.test, 
                                     p_value = out_common_rare$p.value)
        }else{
          resample_common_rare <- Get_Resampling_Pvalue(skat_common_rare)[[1]]
          skat_common_rare <- tibble(SetID_common_rare = SSD.INFO$SetInfo[SET_index,2], 
                                     n_marker = out_common_rare$param$n.marker, 
                                     n_marker_test = out_common_rare$param$n.marker.test, 
                                     emprical_pvalue = resample_common_rare, p_value = out_common_rare$p.value)
        }
        
        return(skat_common_rare)
        
      })
    
    stopCluster(cl)
    
    results <- SKAT_result %>% bind_rows()
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all_common_rare.txt"), col.names = T, row.names = F, sep = "\t")
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
        FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                                  File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE)
        SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
        
        
        
        Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
        out_common_rare <- SKAT_CommonRare(Z, obj, method = "C", CommonRare_Cutoff = 0.05, r.corr.rare = 1,
                                           r.corr.common = 1, test.type = "Common.Only", missing_cutoff = 0.9)
        
        if(re == 0){
          skat_common_rare <- tibble(SetID_common_rare = SSD.INFO$SetInfo[SET_index,2],
                                     n_marker = out_common_rare$param$n.marker,
                                     n_marker_test = out_common_rare$param$n.marker.test, 
                                     p_value = out_common_rare$p.value)
        }else{
          resample_common_rare <- Get_Resampling_Pvalue(skat_common_rare)[[1]]
          skat_common_rare <- tibble(SetID_common_rare = SSD.INFO$SetInfo[SET_index,2], 
                                n_marker = out_common_rare$param$n.marker, 
                                n_marker_test = out_common_rare$param$n.marker.test, 
                                emprical_pvalue = resample_common_rare, p_value = out_common_rare$p.value)
        }
        
        return(skat_common_rare)
        
      })
    
    
    stopCluster(cl)
    
    results <- SKAT_result %>% bind_rows()
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all_NR_gene.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
  }
  
}



### cov making
{
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --indep-pairwise 50 5 0.2 --out skat_pruning", name = "NeuroX"))
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --genome full --out skatQC_IBD", name = "NeuroX")) 
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --read-genome skatQC_IBD.genome --extract skat_pruning.prune.in --mds-plot 20 --cluster --out {name}_MDS",
  #             name = "NeuroX"))
  # 
  # 
  # ###
  # 
  # FAM <- Read_Plink_FAM_Cov(Filename = "../IPDGC.fam", File_Cov = "IPDGC.cov",Is.binary = F)
  # 
  # IPDGC_mds <- fread(file = "IPDGC_MDS.mds", header = T) %>% select(., -SOL)
  # IPDGC_age_control <- fread(file = "IPDGC_control_phenotype.txt", select = c(2,7)) %>% select(., IID = SUBJECT_ID, AGE = Age)
  # IPDGC_age_case <- fread(file = "IPDGC_case_phenotype.txt", select = c(2,7)) %>% select(., IID = SUBJECT_ID, AGE = Age)
  # IPDGC_age <- rbind(IPDGC_age_case, IPDGC_age_control)
  # 
  # IPDGC_cov <- left_join(x = FAM, y = IPDGC_age, by = "IID") %>% left_join(x = ., y = IPDGC_mds, by = c("FID","IID"))
  # 
  # 
  # NeuroX_mds <- fread(file = "NeuroX_MDS.mds", header = T) %>% select(., FID, IID, C1, C2, C3, C4)
  # NeuroX_age <- fread(file = "NeuroX_phenotype.txt", header = T, select = c(2,5)) %>% rename(., IID = SUBJECT_ID) 
  # 
  # IPDGC_cov <- left_join(x = FAM, y = IPDGC_age, by = "IID") %>% left_join(x = ., y = IPDGC_mds, by = c("FID","IID"))
  # NeuroX_cov <- left_join(x = FAM, y = NeuroX_age, by = "IID") %>% left_join(x = ., y = NeuroX_mds, by = c("FID","IID"))
  # 
  # fwrite(x = select(IPDGC_cov, -PID, -MID, -Sex, -Phenotype), file = "IPDGC.cov", row.names = F, sep = "\t")
  # fwrite(x = select(NeuroX_cov, -PID, -MID, -Sex, -Phenotype), file = "NeuroX.cov", row.names = F, sep = "\t")
}  
{
  
  # run_skat_all_resampling <- function(data_name, flag = "default"){
  #   
  #   if(flag == "geneset"){
  #     
  #     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                        File.SetID = paste0(data_name,"_skat.SetID"), 
  #                        File.SSD = paste0(data_name,".SSD"), 
  #                        File.Info = paste0(data_name,".INFO"))
  #     
  #     
  #     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
  #     
  #     if(length(colnames(FAM)) > 20){
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
  #                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #                               C19 + C20, data = FAM, out_type="C", Adjustment = T, n.Resampling = 5000)
  #     }else{
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
  #                             data = FAM, out_type="C", Adjustment = F, n.Resampling = 5000)
  #     }
  #     
  #     # skat run
  #     cl <- makeCluster(detectCores() - 1)
  #     clusterExport(cl, "obj")
  #     # clusterExport(cl, "SSD.INFO")
  #     # clusterEvalQ(cl, obj)
  #     clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
  #     
  #     
  #     SKAT_result <- 
  #       parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index, name = "IPDGC"){
  #         
  #         SSD.INFO <- Open_SSD(File.SSD = paste0(name,".SSD"), File.Info = paste0(name,".INFO"))
  #         Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
  #         out_03 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03)
  #         out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01)
  #         resample_p_03 <- Get_Resampling_Pvalue(out_03)[[1]]
  #         resample_p_01 <- Get_Resampling_Pvalue(out_01)[[1]]
  #         
  #         skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker, 
  #                            n_marker_test = out_03$param$n.marker.test, 
  #                            emprical_pvalue = resample_p_03, p_value = out_03$p.value)
  #         skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker, 
  #                            n_marker_test = out_01$param$n.marker.test, 
  #                            emprical_pvalue = resample_p_01, p_value = out_01$p.value)
  #         
  #         return(bind_cols(skat_001,skat_003))
  #         
  #       })
  #     
  #     stopCluster(cl)
  #     
  #     results <- SKAT_result %>% bind_rows()
  #     
  #     fwrite(x = results, file = paste0(data_name,"_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     Close_SSD()
  #     
  #     
  #     
  #   } else {
  #     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                        File.SetID = paste0(data_name,"_skat_gene.SetID"), 
  #                        File.SSD = paste0(data_name,"_gene.SSD"), 
  #                        File.Info = paste0(data_name,"_gene.INFO"))
  #     
  #     
  #     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
  #     
  #     if(length(colnames(FAM)) > 20){
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
  #                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #                               C19 + C20, data = FAM, out_type="C", Adjustment = F)
  #     }else{
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
  #                             data = FAM, out_type="C", Adjustment = F)
  #     }
  #     
  #     
  #     ## run
  #     out <- list()
  #     out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
  #     out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
  #     
  #     fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     Close_SSD()
  #   }
  #   
  # }
  #   run_skat_gender_age <- function(data_name, flag = "default"){
  #   
  #   if(flag == "geneset"){
  # 
  #       Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                          File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                          File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                          File.SetID = paste0(data_name,"_skat.SetID"), 
  #                          File.SSD = paste0(data_name,".SSD"), 
  #                          File.Info = paste0(data_name,".INFO"))
  #       
  #       
  #       FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                                 File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #       SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
  #       
  #       # if(length(colnames(FAM)) > 20){
  #       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + 
  #       #                           C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #       #                           C19 + C20, data = FAM, out_type="C", Adjustment = F)
  #       # }else{
  #       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4, 
  #       #                         data = FAM, out_type="C", Adjustment = F)
  #       # }
  #       
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE, data = FAM, out_type="C", Adjustment = F)
  #     
  #       ## run
  #       out <- list()
  #       out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
  #       out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
  #         
  #       fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_result.txt"), col.names = T, row.names = F, sep = "\t")
  #       fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_result.txt"), col.names = T, row.names = F, sep = "\t")
  #       Close_SSD()
  #     
  #   } else {
  #       Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                          File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                          File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                          File.SetID = paste0(data_name,"_skat_gene.SetID"), 
  #                          File.SSD = paste0(data_name,"_gene.SSD"), 
  #                          File.Info = paste0(data_name,"_gene.INFO"))
  #       
  #       
  #       FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                                 File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #       SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
  #       
  #       # if(length(colnames(FAM)) > 20){
  #       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + 
  #       #                           C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #       #                           C19 + C20, data = FAM, out_type="C", Adjustment = F)
  #       # }else{
  #       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4, 
  #       #                         data = FAM, out_type="C", Adjustment = F)
  #       # }
  #       
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE, data = FAM, out_type="C", Adjustment = F)
  #       
  #       ## run
  #       out <- list()
  #       out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
  #       out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
  #       
  #       fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result.txt"), col.names = T, row.names = F, sep = "\t")
  #       fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result.txt"), col.names = T, row.names = F, sep = "\t")
  #       Close_SSD()
  #       
  #     }
  # 
  # } # flag (geneset, gene)
  #   run_skat_all_SSD_ALL <- function(data_name, flag = "default"){
  #   
  #   if(flag == "geneset"){
  #     
  #     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                        File.SetID = paste0(data_name,"_skat.SetID"), 
  #                        File.SSD = paste0(data_name,".SSD"), 
  #                        File.Info = paste0(data_name,".INFO"))
  #     
  #     
  #     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
  #     
  #     if(length(colnames(FAM)) > 20){
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
  #                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #                               C19 + C20, data = FAM, out_type="C", Adjustment = F)
  #     }else{
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
  #                             data = FAM, out_type="C", Adjustment = F)
  #     }
  #     
  #     
  #     ## run
  #     out <- list()
  #     out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
  #     out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
  #     
  #     fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     Close_SSD()
  #     
  #   } else {
  #     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
  #                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
  #                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
  #                        File.SetID = paste0(data_name,"_skat_gene.SetID"), 
  #                        File.SSD = paste0(data_name,"_gene.SSD"), 
  #                        File.Info = paste0(data_name,"_gene.INFO"))
  #     
  #     
  #     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
  #                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
  #     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
  #     
  #     if(length(colnames(FAM)) > 20){
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
  #                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
  #                               C19 + C20, data = FAM, out_type="C", Adjustment = F)
  #     }else{
  #       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
  #                             data = FAM, out_type="C", Adjustment = F)
  #     }
  #     
  #     
  #     ## run
  #     out <- list()
  #     out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
  #     out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
  #     
  #     fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
  #     Close_SSD()
  #   }
  #   
  # }
}
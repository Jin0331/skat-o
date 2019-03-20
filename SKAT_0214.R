#library and path
{
  library(dplyr);library(stringr);library(data.table);library(SKAT);library(glue)
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/tool/:/home/jinoo/annovar/:/home/lee/epacts_0913/bin/:/home/jinoo/ensembl-vep")
}

data_name <- c("IPDGC", "NeuroX", "PPMI")

# gzip vcf
{
  
  for(index in 2:3)
    system(glue("gzip -dk {vcf}", vcf = paste0(data_name[index], "_QC.vcf.gz")))
}


setwd("/home/jinoo/skat-o/0311_test/")
for(index in 1:3){
  # vcf to plink format
  print("skatB")
  {
      # system(glue("mkdir {name}", name = data_name[index]))
      setwd(paste0(data_name[index], "/"))
      system(glue("plink --vcf /home/jinoo/skat-o/SKAT_data/{vcf} --recode --double-id --keep-allele-order --out skatB", vcf = paste0(data_name[index],"_QC.vcf")))
      system("plink --file skatB --make-bed --set-missing-var-ids @:# --out skatB")
    
  }
  
  ## vcf_fix load
  print("fix_load")
  {
    test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name[index], "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
      as_tibble()
    test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
    test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
    test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
  }
  
  print("fix_convert")
  ## blank id --> chr:pos
  {
    test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
    for(change in 1:length(test_id)){
      if(test_id[change] == ""){
        test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
      }
      if(change %% 10000 == 0)
        print(change)
    }
    test_fix$ID <- test_id
    
    fwrite(x = data.frame(x=test_fix$ID,y=test_fix$ALT), file = "a1.txt",col.names = F, row.names = F, sep = "\t")
    system("plink --bfile skatB --make-bed --a1-allele a1.txt 2 1 '#' --out skatB_com")
    
    rm(test_ch, test_id, test_pos, change)
  }
  
  #### geneset
  #geneset load
  print("geneset load")
  {
    geneset_merge <- fread(file = "/home/jinoo/skat-o/SKAT_data/1 gene sets 190308 for SKAT_O.txt", header = T, 
                           sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
    
    col_name <- c("GENE","O1","O2","O2M1","O2M2","O2M3","O3","O4","O5","O6", "O7", "O8", "O9","M2","M2P")
    colnames(geneset_merge) <- col_name
    
    one_hot <- lapply(X=select(geneset_merge,-GENE, -M2), FUN = function(x){ifelse(nchar(x) > 0, 1, 0)}) %>% bind_cols()
    geneset_onehot <- cbind(GENE=geneset_merge$GENE, one_hot) %>% mutate(GENE = as.character(GENE)) %>% cbind(M2 = geneset_merge$M2)
    
    for(geneset in 2:ncol(geneset_onehot)){
      for(k in 1:nrow(geneset_onehot)){
        if(geneset_onehot[k, geneset] >= 1){
          geneset_onehot[k,geneset] = geneset_onehot[k,1]
        } else{
           geneset_onehot[k,geneset] = NA
        }
        
      }
    }
    geneset_merge <- geneset_onehot %>% as.list()
    col_name <- names(geneset_merge)
  }
  # making setid_ geneset
  print("gneneset skat.id")
  {
  for(gene in 2:length(col_name)){ #length(col_name)
    geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
    print(col_name[gene])
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
    print(nrow(nonsynonymous))
    ### 2. CADD > 12.37 variant
    cadd <- subset(variant, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                     | ExonicFunc.knownGene == "stopgain"
                                                     | ExonicFunc.knownGene ==  "stoploss"
                                                     | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                     | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                     | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                     | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
    cadd <- subset(cadd, select = "ID")[,1]
    print(nrow(cadd))
    ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
    lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                                    | ExonicFunc.knownGene ==  "stoploss"
                                                    | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                    | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                    | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                    | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
  
    lof <- subset(lof, select = "ID")[,1]
    print(nrow(lof))
    type <- list(nonsynonymous = nonsynonymous, cadd = cadd, lof = lof)
  
    setID <- list()
    for(j in 1:length(type)){
      temp <- data.frame(TYPE=rep(names(type[j]), length(type[[j]])), stringsAsFactors = F)
      setID[[j]] <- cbind(temp, ID=type[[j]])
    }
    setID <- bind_rows(setID)
    setID$TYPE <- paste0(setID$TYPE, "__",col_name[gene])
    # system("rm -rf skat.SetID")
    fwrite(x = setID, "skat.SetID", sep = "\t",row.names = F, quote = F, col.names = F,append = T)
  }
  }
  
  print('geneset skat_gene.id')
  ##gene setID
  {
    variant_gene <- unique(geneset_merge[[1]])
    temp <- list()
    for(i in 1:length(variant_gene)){
      t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED > 12.37))
      temp[[i]] <- subset.data.frame(t1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" 
                                                   # | ExonicFunc.knownGene == "stopgain" | ExonicFunc.knownGene ==  "stoploss" 
                                                   # | ExonicFunc.knownGene ==  "frameshift_deletion"| ExonicFunc.knownGene ==  "frameshift_insertion"
                                                   # | ExonicFunc.knownGene ==  "frameshift_block_substitution"| ExonicFunc.knownGene ==  "splicing"
                                     ), select = c("Gene.knownGene","ID","ID2", "CADD13_PHRED"))
    }
    setID <- bind_rows(temp)
    
    write.table(x = select(setID, Gene.knownGene, ID), "skat_gene.SetID", sep = "\t",row.names = F, quote = F, col.names = F,append = T)
    write.table(x = select(setID, Gene.knownGene, ID2, CADD13_PHRED), "skat_gene_variant.txt", sep = "\t",row.names = F, quote = F, col.names = F,append = T)
  }

  ## pre-processing for skat
  {
    Generate_SSD_SetID(File.Bed = "skatB_com.bed",
                       File.Bim = "skatB_com.bim", 
                       File.Fam = "skatB_com.fam",
                       File.SetID = "skat.SetID", 
                       File.SSD = "skatB_com.SSD", 
                       File.Info = "skatB_com.INFO")
    
    system(glue("cp /home/jinoo/skat-o/SKAT_data/{ped} skatB_com.fam", ped = paste0(data_name[index], ".ped"))) 
    
    FAM<-Read_Plink_FAM(Filename = "skatB_com.fam", Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = "skatB_com.SSD", File.Info = "skatB_com.INFO")
    # Close_SSD()
    obj <-SKAT_Null_Model(Phenotype ~ Sex, data = FAM, out_type="C",Adjustment = F)
    # out <- SKAT.SSD.All(SSD.INFO, obj, method = "SKATO")
    # out <- SKAT(Z, obj, method = "optimal", max_maf = 0.03)
  }

  ## run
  system.time({
    out <- list()
    out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
    out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01)
    
    # system("mkdir result");setwd("result/");system("cp ../skat_gene_variant.txt skat_gene_variant.txt")
    system(glue("mkdir {name}_result", name = data_name[index]));setwd(paste0(data_name[index],"_result/"))
    fwrite(x = out[[1]]$results, file = paste0(data_name[index],"_cut_009.txt"), col.names = T, row.names = F, sep = "\t")
    fwrite(x = out[[2]]$results, file = paste0(data_name[index],"_default.txt"), col.names = T, row.names = F, sep = "\t")
  })
  setwd("/home/jinoo/skat-o/0311_test/")
  gc()
  print(paste0(data_name[index], "_done!!!!!!!!!"))
  
}

## adjusted p-value
{
# 
# for(i in 1:length(file)){
#   temp <- read.table(file = file[i], header = F)
#   system(glue("rm -rf {remove}", remove = file[i]))
#   temp_colname <- c("#CHROM","BEGIN","END","MARKER_ID","NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE","STATRHO");colnames(temp) <- temp_colname
#   PVALUE_ADJUSTED <- p.adjust(temp$PVALUE, p.adjust.methods[4])
#   temp <- cbind(temp, PVALUE_ADJUSTED)
#   write.table(x = temp, file = file[i], col.names = T, row.names = F, quote = F, sep = "\t")
# }
}

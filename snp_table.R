# library path
library_load <- function(){
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC)
  library(tidyverse);library(progress);  library(parallel)
  # tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/tool/:")
}

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
  
  DEG1 <- DEG$DEG1 %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG1) <- "DEG1"
  DEG2 <- DEG$DEG1YJK[DEG$DEG1YJK != ""] %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG2) <- "DEG1YJK"
  DEG3 <- DEG$DEG1Mito[DEG$DEG1Mito != ""] %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG3) <- "DEG1Mito"
  
  geneset_onehot <- geneset_onehot %>% bind_cols(., DEG1) %>% bind_cols(., DEG2) %>% bind_cols(., DEG3)
  
  
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
    mutate(., CHROMPOS = paste(CHROM, POS, sep = ":"))
  
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
  select(clinvar[which(str_detect(clinvar$CHROMPOS, CHROMPOS))[1], ], CLNSIG)[,1] %>%
    return()
}
re_CLNDISDB <- function(CHROMPOS, clinvar){
  select(clinvar[which(str_detect(clinvar$CHROMPOS, CHROMPOS))[1], ], CLNDISDB)[,1] %>%
    return()
}


snp_table <- function(data_name, index_, clinvar){
  core <- detectCores() - 1
  geneset_result <- list()
  geneset <- geneset_load()
  fix <- fix_load(data_name)
  # clinvar <- clinvar_load()
  
  freq <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,"_freq.frq"), header = T) %>%
    rename(ID = SNP)
  fix <- left_join(x = fix, y = freq, by = "ID")
  
  geneset_extract(geneset[[1]], geneset[[2]], fix, data_name, index_)
  
  if(data_name == "IPDGC"){
    for(index in index_){
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_non-syn.txt --recode A --out test_dosage_nonsyn", 
                  data_name = data_name, geneset = geneset[[2]][index]))
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_lof.txt --recode A --out test_dosage_lof", 
                  data_name = data_name, geneset = geneset[[2]][index]))
      
      nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
      lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
      
      print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
      nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
        anno_data_nonsyn <- filter(fix, ID == str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]) %>%
          select(., CHROM, POS, ID, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
          mutate(VUS = re_CLNSIG(paste(CHROM, POS, sep = ":"), clinvar), 
                 CLNDISDB = re_CLNDISDB(paste(CHROM, POS, sep = ":"), clinvar))
        sample_nonsyn <- select(nonsyn_dosage, c(1:2, col_len)) %>% filter(., .[,3] >= 1) %>% select(., IID, PHENOTYPE)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_nonsyn) == 0) return(anno_data_mul)
        
        for(i in 1:nrow(sample_nonsyn))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_nonsyn)  
        return(bind_cols(sample_nonsyn, anno_data_mul))
      }, mc.cores = core) %>% bind_rows()
      
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
        anno_data_lof <- filter(fix, ID == str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]) %>%
          select(., CHROM, POS, ID, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
          mutate(VUS = re_CLNSIG(paste(CHROM, POS, sep = ":"), clinvar), 
                 CLNDISDB = re_CLNDISDB(paste(CHROM, POS, sep = ":"), clinvar))
        sample_lof <- select(lof_dosage, c(1:2, col_len)) %>% filter(., .[,3] >= 1) %>% select(., IID, PHENOTYPE)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_lof) == 0) return(anno_data_mul)
        
        for(i in 1:nrow(sample_lof))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_lof)    
        return(bind_cols(sample_lof, anno_data_mul))
      }, mc.cores = core) %>% bind_rows()
      
      nonsyn_result <- nonsyn_result %>%
        mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
               Datasets = data_name, geneset = geneset[[2]][index], 
               Func = "Nonsynonymous")
      
      lof_result <- lof_result %>%
        mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
               Datasets = data_name, 
               geneset = geneset[[2]][index], Func = "LoF")
      
      geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
    }
    
  }else{ # NeuroX
    for(index in index_){
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_non-syn.txt --recode A --out test_dosage_nonsyn", 
                  data_name = data_name, geneset = geneset[[2]][index]))
      system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{data_name} --extract {data_name}_{geneset}_lof.txt --recode A --out test_dosage_lof", 
                  data_name = data_name, geneset = geneset[[2]][index]))
      
      nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)
      lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX)

      print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
      nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
        anno_data_nonsyn <- filter(fix, ID == str_sub(colnames(nonsyn_dosage)[col_len], end = -3)) %>%
          select(., CHROM, POS, ID, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
          mutate(VUS = re_CLNSIG(paste(CHROM, POS, sep = ":"), clinvar), 
                 CLNDISDB = re_CLNDISDB(paste(CHROM, POS, sep = ":"), clinvar))
        sample_nonsyn <- select(nonsyn_dosage, c(1:2, col_len)) %>% filter(., .[,3] >= 1) %>% select(., IID, PHENOTYPE)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_nonsyn) == 0) return(anno_data_mul)
        
        for(i in 1:nrow(sample_nonsyn))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_nonsyn)  
        return(bind_cols(sample_nonsyn, anno_data_mul))
      }, mc.cores = core) %>% bind_rows()

      
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
        anno_data_lof <- filter(fix, ID == str_sub(colnames(lof_dosage)[col_len], end = -3)) %>%
          select(., CHROM, POS, ID, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
          mutate(VUS = re_CLNSIG(paste(CHROM, POS, sep = ":"), clinvar), 
                 CLNDISDB = re_CLNDISDB(paste(CHROM, POS, sep = ":"), clinvar))
        sample_lof <- select(lof_dosage, c(1:2, col_len)) %>% filter(., .[,3] >= 1) %>% select(., IID, PHENOTYPE)
        
        anno_data_mul <- tibble(.rows = 0)
        if(nrow(sample_lof) == 0) return(anno_data_mul)
        
        for(i in 1:nrow(sample_lof))
          anno_data_mul <- bind_rows(anno_data_mul, anno_data_lof)  
        return(bind_cols(sample_lof, anno_data_mul) )
      }, mc.cores = core) %>% bind_rows()
      
      nonsyn_result <- nonsyn_result %>%
        mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
               Datasets = data_name, geneset = geneset[[2]][index], 
               Func = "Nonsynonymous")
      
      lof_result <- lof_result %>%
        mutate(., AnnotationClinvar = ifelse(paste(CHROM,POS,sep = ":") %in% clinvar$CHROMPOS, "TRUE","FALSE"),
               Datasets = data_name, 
               geneset = geneset[[2]][index], 
               Func = "LoF")
      
      geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
    }
  }
  
  return(geneset_result)
}



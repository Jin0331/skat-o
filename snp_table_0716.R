# library path
library_load <- function(){
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(tidyverse);library(progress);  library(parallel)
  library(tidyselect);library(magrittr)
  
  # tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/tool/:")
}

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
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
    
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
    
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
    
    
    return(test_fix)
    
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
    
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
    test_fix$Gene.knownGene <- str_replace(test_fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
    
    
    return(test_fix)
  }
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
    
    nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX) %>% as_tibble()
    lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX) %>% as_tibble()
    
    nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
    # nonsyn_dosage[which(nonsyn_dosage$IID == "UMARY-1571"), which(colnames(nonsyn_dosage) == "rs138008832_A")] <- 1 ######### remove
    
    lof_dosage[is.na(lof_dosage)] <- 0
    
    print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
    nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
      anno_data_nonsyn <- filter(fix, ID == str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]) %>%
        select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
        mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
      
      sample_nonsyn <- select(nonsyn_dosage, c(1:2, col_len)) %>% 
        # filter(., .[,3] >= 1) %>% 
        mutate(Heterozygous_M = ifelse(.[,3] == 1, TRUE, FALSE),
               Homozygous_M = ifelse(.[,3] == 2, TRUE, FALSE),
               Heterozygous = ifelse(.[,3] == 0, TRUE, FALSE)) %>%
        select(., IID, PHENOTYPE, Heterozygous, Heterozygous_M, Homozygous_M)
      
      anno_data_mul <- tibble(.rows = 0)
      if(nrow(sample_nonsyn) == 0) return(anno_data_mul)
      
      for(i in 1:nrow(sample_nonsyn))
        anno_data_mul <- bind_rows(anno_data_mul, anno_data_nonsyn)  
      return(bind_cols(sample_nonsyn, anno_data_mul))
    }, mc.cores = core) %>% bind_rows()
    
    print(paste(geneset[[2]][index],"lof", sep = " "))
    lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
      anno_data_lof <- filter(fix, ID == str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]) %>%
        select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>% 
        mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
      
      sample_lof <- select(lof_dosage, c(1:2, col_len)) %>% 
        # filter(., .[,3] >= 1) %>% 
        mutate(Heterozygous_M = ifelse(.[,3] == 1, TRUE, FALSE),
               Homozygous_M = ifelse(.[,3] == 2, TRUE, FALSE),
               Heterozygous = ifelse(.[,3] == 0, TRUE, FALSE)) %>%
        select(., IID, PHENOTYPE, Heterozygous, Heterozygous_M, Homozygous_M)
      
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
    
    geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
  }
  
  system("rm -rf *.log");system("rm -rf *.raw");system("rm -rf *.hh")
  system("rm -rf *.txt")
  return(geneset_result)
}
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
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>%
        bind_rows()
      
      
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
        anno_data_lof <- filter(fix, ID == str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]) %>%
          select(., CHROM, POS, ID, REF, ALT, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))}, mc.cores = core) %>% 
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
          mutate(Clinical_Significance = re_CLNSIG(paste(CHROM, ":", POS, REF, ALT), clinvar))},mc.cores = core) %>%
        bind_rows()
      
      
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
        anno_data_lof <- filter(fix, ID == str_sub(colnames(lof_dosage)[col_len], end = -3)) %>%
          select(., CHROM, POS, REF, ALT, ID, Gene.knownGene, AAChange.knownGene, CADD13_PHRED, MAF) %>%
          mutate(Clinical_Significance = re_CLNSIG(paste(CHROM, ":", POS, REF, ALT), clinvar))}, mc.cores = core) %>% 
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
    result[1,1] <- 445 - result[,1] %>% sum
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
table3_ver_3_sample <- function(data_list, geneset_name){
  index <- which(names(data_list) == geneset_name)
  temp <- data_list[[index]]
  temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                            pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                            replacement = "Pathogenic_ALL")
  temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
  sample <- temp$IID %>% unique() %>% sort()
  
  iid <- sample[which(sample == "UMARY-1903")]
  
  sample_result <- mclapply(X = sample, FUN = function(iid){
    # result <- list()
    PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
      unique() 
    PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
    
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
    
    
    tibble(Subject_ID = iid, PHENOTYPE, Pathogenic_ALL = (Pathogenic_Hetro + (Pathogenic_Homo * 2)),
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



{
  # table3_ver_1 <- function(data_list, type = "Pathogenic", geneset_name, MAF_value = 0.03, CADD_score = 20){
  #   index <- which(names(data_list) == geneset_name)
  #   
  #   temp <- data_list[[index]]
  #   temp$Clinical_Significance <- str_replace(temp$Clinical_Significance,
  #                                             pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)",
  #                                             replacement = "Pathogenic_ALL")
  #   temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
  #   
  #   if(type == "Pathogenic"){
  #     temp2 <- temp %>% arrange(IID) %>% filter(MAF <= MAF_value, Clinical_Significance == "Pathogenic_ALL") %>%
  #       table3_CLSIG() %>% bind_rows()
  #     temp2[is.na(temp2)] <- 0
  #     
  #     control <- temp2 %>% filter(PHENOTYPE == 1) %>%
  #       rename(Pathogenic_ALL_control = Pathogenic_ALL) %>%
  #       table3_calc() 
  #     case <- temp2 %>% filter(PHENOTYPE == 2) %>%
  #       rename(Pathogenic_ALL_case = Pathogenic_ALL) %>%
  #       table3_calc() 
  #     
  #   } else if(type == "LoF"){
  #     temp2 <- temp %>% arrange(IID) %>% filter(MAF <= MAF_value, Func == "LoF", Clinical_Significance != "Pathogenic_ALL") %>%
  #       table3_CLSIG() %>% bind_rows()
  #     temp2[is.na(temp2)] <- 0
  #     
  #     lof_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(lof_sum) <- "LoF"
  #     temp2 <- bind_cols(temp2[, 1:2], lof_sum)
  #     
  #     control <- temp2 %>% filter(PHENOTYPE == 1) %>%
  #       rename(LoF_control = LoF) %>%
  #       table3_calc() 
  #     case <- temp2 %>% filter(PHENOTYPE == 2) %>%
  #       rename(LoF_case = LoF) %>%
  #       table3_calc() 
  #     
  #   }else if(type == "CADD_MAF"){
  #     temp2 <- temp %>% arrange(IID) %>% filter(MAF <= MAF_value, CADD13_PHRED <= CADD_score, Clinical_Significance != "Pathogenic_ALL") %>%
  #       table3_CLSIG() %>% bind_rows()
  #     temp2[is.na(temp2)] <- 0
  #     
  #     cadd_maf_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2]
  #     colnames(cadd_maf_sum) <- "cadd_maf"
  #     temp2 <- bind_cols(temp2[, 1:2], cadd_maf_sum)
  #     
  #     control <- temp2 %>% filter(PHENOTYPE == 1) %>%
  #       rename(CADD_MAF_control = cadd_maf) %>%
  #       table3_calc() 
  #     
  #     case <- temp2 %>% filter(PHENOTYPE == 2) %>%
  #       rename(CADD_MAF_case = cadd_maf) %>%
  #       table3_calc() 
  #     
  #     
  #     
  #   } else if(type == "ALL"){ # not provide
  #     temp2 <- temp %>% arrange(IID) %>% filter(MAF <= 0.03) %>%
  #       table3_CLSIG() %>% bind_rows()
  #     temp2[is.na(temp2)] <- 0
  #     
  #     all_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(all_sum) <- "ALL"
  #     temp2 <- bind_cols(temp2[, 1:2], all_sum)
  #     
  #     control <- temp2 %>% filter(PHENOTYPE == 1) %>%
  #       rename(ALL_control = ALL) %>%
  #       table3_calc() 
  #     case <- temp2 %>% filter(PHENOTYPE == 2) %>%
  #       rename(ALL_case = ALL) %>%
  #       table3_calc() 
  #   }
  #   
  #   
  #   result <- bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), control, case) %>% 
  #     column_to_rownames(., var = "rowname")
  #   
  #   if(result[1,1] == 0)
  #     result[1,1] <- 333 - result[,1] %>% sum
  #   if(result[1,2] == 0)
  #     result[1,2] <- 445 - result[,2] %>% sum
  #   
  #   return(result)
  # }
  
  
  # re_CLNDISDB <- function(CHROMPOS, clinvar){
  #   select(clinvar[which(str_detect(clinvar$CHROMPOS, CHROMPOS))[1], ], CLNDISDB)[,1] %>%
  #     return()
  # }
  
  #else if(type == "ALL"){ # not provide
  #   temp2 <- temp %>% arrange(IID) %>% filter(MAF <= 0.03, Heterozygous != TRUE) %>%
  #     table3_CLSIG() %>% bind_rows()
  #   temp2[is.na(temp2)] <- 0
  # 
  #   all_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(all_sum) <- "ALL"
  #   temp2 <- bind_cols(temp2[, 1:2], all_sum)
  # 
  #   control <- temp2 %>% filter(PHENOTYPE == 1) %>%
  #     rename(ALL_control = ALL) %>%
  #     table3_calc()
  #   case <- temp2 %>% filter(PHENOTYPE == 2) %>%
  #     rename(ALL_case = ALL) %>%
  #     table3_calc()
  # }
  
  
  # for(iid in sample){
  #   
  #   PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
  #     unique() 
  #   PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
  #   
  #   # Pathogenic
  #   {
  #     Pathogenic_Hetro <- temp %>% 
  #       filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
  #              IID == iid, Heterozygous_M == T) %>% 
  #       select(Clinical_Significance) %>%
  #       group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
  #     
  #     if(Pathogenic_Hetro != 0){
  #       Hetero_name <- temp %>% 
  #         filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
  #                IID == iid, Heterozygous_M == T) %>% .$ID
  #       Pathogenic_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
  #     } else {Pathogenic_Hetro_ID <- " "}
  #     
  #     
  #     Pathogenic_Homo <- temp %>% 
  #       filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
  #              IID == iid, Homozygous_M == T) %>% 
  #       select(Clinical_Significance) %>%
  #       group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
  #     
  #     if(Pathogenic_Homo != 0){
  #       Homo_name <- temp %>% 
  #         filter(MAF <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
  #                IID == iid, Homozygous_M == T) %>% .$ID
  #       Pathogenic_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
  #     } else {Pathogenic_Homo_ID <- " "}
  #   }
  #   
  #   # LoF
  #   {
  #     LoF_Hetero <- temp %>% 
  #       filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
  #              IID == iid, Func == "LoF", Heterozygous_M == T) %>%
  #       select(Clinical_Significance) %>%
  #       group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
  #     
  #     if(LoF_Hetero != 0){
  #       Hetero_name <- temp %>% 
  #         filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
  #                IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$ID
  #       LoF_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
  #     } else {LoF_Hetro_ID <- " "}
  #     
  #     LoF_Homo <- temp %>% 
  #       filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
  #              IID == iid, Func == "LoF", Homozygous_M == T) %>%
  #       select(Clinical_Significance) %>%
  #       group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
  #     
  #     if(LoF_Homo != 0){
  #       Homo_name <- temp %>% 
  #         filter(MAF <= 0.03, Clinical_Significance != "Pathogenic_ALL", 
  #                IID == iid, Func == "LoF", Homozygous_M == T) %>% .$ID
  #       LoF_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
  #     } else {LoF_Homo_ID <- " "}
  #   }
  #   
  #   # CADD
  #   {
  #     for(value in c(0.03, 0.01)){
  #       CADD_Hetero <- temp %>% 
  #         filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
  #                IID == iid, Heterozygous_M == T) %>%
  #         select(Clinical_Significance) %>%
  #         group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
  #       
  #       if(CADD_Hetero != 0){
  #         Hetero_name <- temp %>% 
  #           filter(MAF <= value, CADD13_PHRED >= CADD_score, Clinical_Significance != "Pathogenic_ALL", 
  #                  IID == iid, Heterozygous_M == T) %>% .$ID
  #         CADD_Hetro_name <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
  #       } else {CADD_Hetro_name <- " "}
  #       
  #       CADD_Homo <- temp %>% 
  #         filter(MAF <= value, Clinical_Significance != "Pathogenic_ALL", 
  #                IID == iid, Homozygous_M == T) %>%
  #         select(Clinical_Significance) %>%
  #         group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
  #       
  #       if(CADD_Homo != 0){
  #         Homo_name <- temp %>% 
  #           filter(MAF <= value, Clinical_Significance != "Pathogenic_ALL", 
  #                  IID == iid, Homozygous_M == T) %>% .$ID
  #         CADD_Homo_name <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
  #       } else {CADD_Homo_name <- " "}
  #       
  #       if(value == 0.03){
  #         CADD_MAF003_Hetero <- CADD_Hetero
  #         CADD_MAF003_Hetero_ID <- CADD_Hetro_name
  #         CADD_MAF003_Homo <- CADD_Homo
  #         CADD_MAF003_Homo_ID <- CADD_Homo_name
  #       } else{
  #         CADD_MAF001_Hetero <- CADD_Hetero
  #         CADD_MAF001_Hetero_ID <- CADD_Hetro_name
  #         CADD_MAF001_Homo <- CADD_Homo
  #         CADD_MAF001_Homo_ID <- CADD_Homo_name
  #       }
  #       
  #     }
  #   }
  #   
  #   
  #   result[[`iid`]] <- tibble(Subject_ID = iid, PHENOTYPE, Pathogenic_ALL = (Pathogenic_Hetro + (Pathogenic_Homo * 2)),
  #                             LoF_ALL = (LoF_Hetero + (LoF_Homo * 2)),
  #                             CADD_MAF003_ALL = (CADD_MAF003_Hetero + (CADD_MAF003_Homo * 2)),
  #                             CADD_MAF001_ALL = (CADD_MAF001_Hetero + (CADD_MAF001_Homo * 2)),
  #                             Pathogenic_Hetro, Pathogenic_Hetro_ID,
  #                             Pathogenic_Homo, Pathogenic_Homo_ID,
  #                             LoF_Hetero, LoF_Hetro_ID, LoF_Homo, LoF_Homo_ID,
  #                             CADD_MAF003_Hetero, CADD_MAF003_Hetero_ID,
  #                             CADD_MAF003_Homo, CADD_MAF003_Homo_ID,
  #                             CADD_MAF001_Hetero, CADD_MAF001_Hetero_ID,
  #                             CADD_MAF001_Homo, CADD_MAF001_Homo_ID) 
  #   
  #   print(iid) ## remove
  #   
  # }
  
}

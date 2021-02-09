library_load <- function(){
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(tidyverse);library(parallel)
  library(tidyselect);library(magrittr);library(SKAT);library(progress);library(RMySQL)
  #library(MetaSKAT)
  # tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/dblab/tool/:")
  
  # power
  data("SKAT.haplotypes")
  attach(SKAT.haplotypes)
}
geneset_KFPD<- function(path){
  print("Geneset load GO Client geneset!")
  # 92 gene_SLC25A4 remove
  geneset_merge <- fread(path, header = T, stringsAsFactors = F, data.table = F) %>% as_tibble()
  return_list <- list()
  return_list[[1]] <- lapply(X = geneset_merge, FUN = function(value){
    value %>% .[. != ""] %>% return()
  })
  return_list[[2]] <- colnames(geneset_merge)
  
  return_list[[3]] <- lapply(X = return_list[[1]], function(x){length(x)}) %>% unlist()
  return(return_list)
} 
fix_load_KFPD <- function(path, plink_path){
  print("fix load for Korean FPD")
  fix <- fread(path) %>% 
    as_tibble() 
  fix$Gene.knownGene <- str_replace(fix$Gene.knownGene, pattern = "^PARK2$", replacement = "PRKN")
  fix$Gene.knownGene <- str_replace(fix$Gene.knownGene, pattern = "^C10orf2$", replacement = "TWNK")
  
  if("ID2" %in% colnames(fix) == F){
    id2 <- fread(paste(plink_path, "bim", sep = ".")) %>% as_tibble() %>% select(ID2 = V2)
    fix <- fix %>% bind_cols(id2, .)
  }
  
  return(fix)  
}

clinvar_load <- function(clinvar_path){
  # clinvar preprocessing
  clinvar <- fread(file = clinvar_path, header = T) %>%  rename(CHROM = `#CHROM`) %>% 
    mutate(., CHROM = str_replace(string = CHROM, pattern = "X", replacement = "23")) %>%
    mutate(., CHROM = str_replace(string = CHROM, pattern = "Y", replacement = "24")) %>%
    mutate(., CHROM = str_replace(string = CHROM, pattern = "MT", replacement = "26")) %>%
    mutate(., CHROMPOS = paste(CHROM, POS, sep = ":")) %>%
    mutate(., CHROMPOS_R_A = paste0(CHROMPOS, REF, ALT))
  
  info <- mclapply(X = clinvar$INFO, FUN = clinvar_info, mc.cores = detectCores() - 1) %>% bind_rows()
  clinvar %>% select(-INFO) %>% mutate(CHROM = as.integer(CHROM)) %>%  bind_cols(info) %>% 
    as_tibble() %>% return()
}
clinvar_info <- function(temp){
  str_split(temp,pattern = ";")[[1]] %>% 
    as_tibble() %>% separate(value, c("info","value"), sep = "=") %>% 
    spread(key = "info", value = "value") %>% return()
}

re_CLNSIG <- function(CHROMPOS, clinvar){
  select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], CLNSIG) %>% pull(1) %>%
    return()
}
re_CLNSIGCONF <- function(CHROMPOS, clinvar){
  select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], CLNSIGCONF) %>% pull(1) %>%
    return()
}
re_CLINVAR_ID <- function(CHROMPOS, clinvar){
  select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], ID) %>% pull(1) %>%
    return()
}
re_CLNHGVS <- function(CHROMPOS, clinvar){
  select(clinvar[which(str_detect(clinvar$CHROMPOS_R_A, paste0("^",CHROMPOS, "$")))[1], ], CLNHGVS) %>% pull(1) %>%
    return()
}
snp_table_func1 <- function(table){
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
snp_table_func2 <- function(table_list){
  result_list <- list()
  col_name <- colnames(table_list)
  
  # if(nrow(table_list) == 0) return(tibble())
  
  for(index in 3:ncol(table_list)){
    temp <- table_list %>% select(col_name[index]) %>% 
      group_by_all() %>% 
      summarise(value = n(), .groups = "keep") %>%
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

geneset_extract <- function(geneset_merge, col_name, fix, data_name, index_){
  print("geneset SetID making!!")
  
  for(gene in index_){ # 2:length(col_name), c(2,3,7,16,17)
    geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
    variant <- list()
    ## Parkinson geneset variant
    for(i in 1:length(geneset)){
      variant[[i]] <- subset(fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")),
                             select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","PHRED","ID2"))
    }
    variant <- bind_rows(variant)
    
    ### 1. nonsynonymous geneset
    nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"))
    fwrite(x = subset(nonsynonymous, select = "ID2")[,1], file = paste0(data_name,"_",col_name[gene],"_non-syn.txt"), col.names = F)
    
    
    
    ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
    lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                      | ExonicFunc.knownGene ==  "stoploss"
                                      | ExonicFunc.knownGene ==  "frameshift_deletion"
                                      | ExonicFunc.knownGene ==  "frameshift_insertion"
                                      | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                      | ExonicFunc.knownGene ==  "splicing"))
    
    fwrite(x = subset(lof, select = "ID2")[,1], file = paste0(data_name,"_",col_name[gene],"_lof.txt"), col.names = F)
  }
} 
snp_table_variant <- function(data_name, index_, clinvar, fix, geneset, plink_path){
  core <- detectCores() - 1
  geneset_result <- list()
  
  if(!dir.exists("geneset_extract") | getwd() != "geneset_extract"){
    dir.create("geneset_extract")
  } 
  
  setwd("geneset_extract/")
  geneset_extract(geneset[[1]], geneset[[2]], fix, data_name, index_)
  
  for(index in index_){
    system(glue("/tools/plink --bfile {plink_path} --extract {data_name}_{geneset}_non-syn.txt --allow-no-sex --recode A --out test_dosage_nonsyn",
                data_name = data_name, geneset = geneset[[2]][index], plink_path = plink_path), ignore.stdout = T)
    system(glue("/tools/plink --bfile {plink_path} --extract {data_name}_{geneset}_lof.txt --allow-no-sex --recode A --out test_dosage_lof",
                data_name = data_name, geneset = geneset[[2]][index], plink_path = plink_path), ignore.stdout = T)
    
    nonsyn_dosage <- tibble()
    tryCatch(expr = nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX) %>% as_tibble(),
             error = function(e){ print(e) })
    
    if(nrow(nonsyn_dosage) > 0){
      nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
      print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
      nonsyn_result <- mclapply(X = 3:ncol(nonsyn_dosage), FUN = function(col_len){
        
        if(!str_detect(string = colnames(nonsyn_dosage)[col_len], pattern = "^rs")){
          id_split <- str_split(colnames(nonsyn_dosage)[col_len], "_") %>% unlist()
          temp <- paste(id_split[1], id_split[2], id_split[3], sep = "_")
        } else {
          temp <- str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]
        }
        
        anno_data_nonsyn <- filter(fix, ID2 == temp) %>%
          select(., CHROM, POS, ID2, REF, ALT, FILTER, Gene.knownGene, AAChange.knownGene, PHRED, AF, AF_eas, AF_eas_kor) %>% 
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLINVAR_ID = re_CLINVAR_ID(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLNHGVS = re_CLNHGVS(paste0(CHROM, ":", POS, REF, ALT), clinvar))},
        mc.cores = core) %>% bind_rows()
      
      nonsyn_result <- nonsyn_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               geneset = geneset[[2]][index], Func = "Nonsynonymous")
    } else {
      nonsyn_result <- NULL
    } # nrows of nonsyn_dosage > 0
    
    lof_dosage <- tibble()
    tryCatch(expr = lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX) %>% as_tibble(),
             error = function(e){ print(e) })
    
    if(nrow(lof_dosage) > 0){
      lof_dosage[is.na(lof_dosage)] <- 0
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 3:ncol(lof_dosage), FUN = function(col_len){
        
        if(!str_detect(string = colnames(lof_dosage)[col_len], pattern = "^rs")){
          id_split <- str_split(colnames(lof_dosage)[col_len], "_") %>% unlist()
          temp <- paste(id_split[1], id_split[2], id_split[3], sep = "_")
        } else {
          temp <- str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]
        }
        
        anno_data_lof <- filter(fix, ID2 == temp) %>%
          select(., CHROM, POS, ID2, REF, ALT, FILTER, Gene.knownGene, AAChange.knownGene, PHRED, AF, AF_eas, AF_eas_kor) %>%
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLNSIGCONF = re_CLNSIGCONF(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLINVAR_ID = re_CLINVAR_ID(paste0(CHROM, ":", POS, REF, ALT), clinvar),
                 CLNHGVS = re_CLNHGVS(paste0(CHROM, ":", POS, REF, ALT), clinvar)
          )},mc.cores = core) %>%
        bind_rows()
      
      lof_result <- lof_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT)%in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               geneset = geneset[[2]][index], Func = "LoF")
    } else {
      lof_result <- NULL
    }
    
    geneset_result[[geneset[[2]][index]]] <- bind_rows(nonsyn_result, lof_result)
  }
  
  return(geneset_result)
}
snp_table_preprocessing <- function(data_name, index_, clinvar, fix, geneset, plink_path){
  core <- detectCores() - 1
  geneset_result <- list()
  
  if(!dir.exists("geneset_extract")){
    dir.create("geneset_extract")
  } 
  
  setwd("geneset_extract/")
  geneset_extract(geneset[[1]], geneset[[2]], fix, data_name, index_)
  
  for(index in index_){
    system(glue("/tools/plink --bfile {plink_path} --extract {data_name}_{geneset}_non-syn.txt --allow-no-sex --recode A --out test_dosage_nonsyn",
                data_name = data_name, geneset = geneset[[2]][index], plink_path = plink_path), ignore.stdout = T)
    system(glue("/tools/plink --bfile {plink_path} --extract {data_name}_{geneset}_lof.txt  --allow-no-sex --recode A --out test_dosage_lof",
                data_name = data_name, geneset = geneset[[2]][index], plink_path = plink_path), ignore.stdout = T)
    
    nonsyn_dosage <- tibble()
    tryCatch(expr = nonsyn_dosage <- fread(file = "test_dosage_nonsyn.raw", header = T) %>% select(-FID, -PAT, -MAT) %>% as_tibble(),
             error = function(e){ print(e) })
    
    if(nrow(nonsyn_dosage) > 0){
      nonsyn_dosage[is.na(nonsyn_dosage)] <- 0
      print(paste(geneset[[2]][index],"nonsynonymous", sep = " "))
      nonsyn_result <- mclapply(X = 4:ncol(nonsyn_dosage), FUN = function(col_len){
        
        if(!str_detect(string = colnames(nonsyn_dosage)[col_len], pattern = "^rs")){
          id_split <- str_split(colnames(nonsyn_dosage)[col_len], "_") %>% unlist()
          temp <- paste(id_split[1], id_split[2], id_split[3], sep = "_")
        } else {
          temp <- str_split(colnames(nonsyn_dosage)[col_len], pattern = "_")[[1]][1]
        }
        
        anno_data_nonsyn <- filter(fix, ID2 == temp) %>%
          select(., CHROM, POS, ID2, REF, ALT, FILTER, Gene.knownGene, AAChange.knownGene, PHRED, AF, AF_eas, AF_eas_kor) %>% 
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
        
        sample_nonsyn <- select(nonsyn_dosage, c(1:3, col_len)) %>% 
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
      
      nonsyn_result <- nonsyn_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               Datasets = data_name, geneset = geneset[[2]][index], 
               Func = "Nonsynonymous")
    } else {
      nonsyn_result <- NULL
    }
    
    lof_dosage <- tibble()
    tryCatch(expr = lof_dosage <- fread(file = "test_dosage_lof.raw", header = T) %>% select(-FID, -PAT, -MAT) %>% as_tibble(),
             error = function(e){ print(e) })
    
    if(nrow(lof_dosage) > 0){
      lof_dosage[is.na(lof_dosage)] <- 0
      print(paste(geneset[[2]][index],"lof", sep = " "))
      lof_result <- mclapply(X = 4:ncol(lof_dosage), FUN = function(col_len){
        if(!str_detect(string = colnames(lof_dosage)[col_len], pattern = "^rs")){
          id_split <- str_split(colnames(lof_dosage)[col_len], "_") %>% unlist()
          temp <- paste(id_split[1], id_split[2], id_split[3], sep = "_")
        } else {
          temp <- str_split(colnames(lof_dosage)[col_len], pattern = "_")[[1]][1]
        }
        
        anno_data_lof <- filter(fix, ID2 == temp) %>%
          select(., CHROM, POS, ID2, REF, ALT, FILTER, Gene.knownGene, AAChange.knownGene, PHRED, AF, AF_eas, AF_eas_kor) %>% 
          mutate(Clinical_Significance = re_CLNSIG(paste0(CHROM, ":", POS, REF, ALT), clinvar))
        
        sample_lof <- select(lof_dosage, c(1:3, col_len)) %>% 
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
      
      lof_result <- lof_result %>%
        mutate(., AnnotationClinvar = ifelse(paste0(CHROM, ":", POS, REF, ALT) %in% clinvar$CHROMPOS_R_A, "TRUE","FALSE"),
               Datasets = data_name, 
               geneset = geneset[[2]][index], Func = "LoF")
    } else {
      lof_result <- NULL
    }
    
    
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
  setwd("../")
  
  return(geneset_result)
}

snp_table_sample_count <- function(data_list, type = "Pathogenic", geneset_name, MAF_value = 0.03, CADD_score = 20, data_name,
                                   case_number, control_number){
  index <- which(names(data_list) == geneset_name)
  temp <- data_list[[index]]
  temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                            pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                            replacement = "Pathogenic_ALL")
  temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
  
  if(type == "Pathogenic"){
    temp2 <- temp %>% arrange(IID) %>% 
      filter(AF_eas <= MAF_value, Heterozygous != TRUE, Clinical_Significance == "Pathogenic_ALL") %>%
      snp_table_func1() %>% bind_rows()
    if(nrow(temp2) < 1){
      return(NULL)
    }
    
    temp2[is.na(temp2)] <- 0
    control <- temp2 %>% filter(PHENOTYPE == 1) %>%
      rename(Pathogenic_ALL_control = Pathogenic_ALL) %>%
      snp_table_func2() 
    case <- temp2 %>% filter(PHENOTYPE == 2) %>%
      rename(Pathogenic_ALL_case = Pathogenic_ALL) %>%
      snp_table_func2() 
    
  } else if(type == "LoF"){
    temp2 <- temp %>% arrange(IID) %>% 
      filter(Func == "LoF", Heterozygous != TRUE) %>%
      snp_table_func1() %>% bind_rows()
    if(nrow(temp2) < 1){
      return(NULL)
    }
    
    temp2[is.na(temp2)] <- 0
    lof_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(lof_sum) <- "LoF"
    temp2 <- bind_cols(temp2[, 1:2], lof_sum)
    
    control <- temp2 %>% filter(PHENOTYPE == 1) %>%
      rename(LoF_control = LoF) %>%
      snp_table_func2() 
    case <- temp2 %>% filter(PHENOTYPE == 2) %>%
      rename(LoF_case = LoF) %>%
      snp_table_func2() 
    
  }else if(type == "CADD_MAF"){
    temp2 <- temp %>% arrange(IID) %>% 
      filter(AF_eas <= MAF_value, Heterozygous != TRUE, PHRED >= CADD_score) %>%
      snp_table_func1() %>% bind_rows()
    if(nrow(temp2) < 1){
      return(NULL)
    }
    
    temp2[is.na(temp2)] <- 0
    cadd_maf_sum <- temp2[, 3:ncol(temp2)] %>% apply(., MARGIN = 1, sum) %>% enframe() %>% .[2];colnames(cadd_maf_sum) <- "cadd_maf"
    temp2 <- bind_cols(temp2[, 1:2], cadd_maf_sum)
    
    control <- temp2 %>% filter(PHENOTYPE == 1) %>%
      rename(CADD_MAF_control = cadd_maf) %>%
      snp_table_func2() 
    
    case <- temp2 %>% filter(PHENOTYPE == 2) %>%
      rename(CADD_MAF_case = cadd_maf) %>%
      snp_table_func2() 
  } 
  result <- bind_cols(case, control) 
  
  if(result[1,1] == 0)
    result[1,1] <- case_number - result[,1] %>% sum
  if(result[1,2] == 0)
    result[1,2] <- control_number - result[,2] %>% sum  
  
  
  return(result)
}
snp_table_per_sample <- function(data_list, geneset_name, CADD_score = 20){
  index <- which(names(data_list) == geneset_name)
  temp <- data_list[[index]]
  temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                            pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                            replacement = "Pathogenic_ALL")
  temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
  sample <- temp$IID %>% unique() %>% sort()
  
  sample_result <- mclapply(X = sample, FUN = function(iid){
    PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
      unique() 
    PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
    SEX <- temp %>% filter(IID == iid) %>% select(SEX) %>% pull(1) %>%
      unique()
    
    # Pathogenic
    {
      Pathogenic_Hetro <- temp %>% 
        filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Heterozygous_M == T) %>% 
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
      
      if(Pathogenic_Hetro != 0){
        Hetero_name <- temp %>% 
          filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Heterozygous_M == T) %>% .$ID
        Pathogenic_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
      } else {Pathogenic_Hetro_ID <- " "}
      
      
      Pathogenic_Homo <- temp %>% 
        filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Homozygous_M == T) %>% 
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
      
      if(Pathogenic_Homo != 0){
        Homo_name <- temp %>% 
          filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Homozygous_M == T) %>% .$ID
        Pathogenic_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
      } else {Pathogenic_Homo_ID <- " "}
      }
    
    # LoF
    {
      LoF_Hetero <- temp %>% 
        filter(IID == iid, Func == "LoF", Heterozygous_M == T) %>%
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      if(LoF_Hetero != 0){
        Hetero_name <- temp %>% 
          filter(IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$ID
        LoF_Hetro_ID <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
      } else {LoF_Hetro_ID <- " "}
      
      LoF_Homo <- temp %>% 
        filter(IID == iid, Func == "LoF", Homozygous_M == T) %>%
        select(Clinical_Significance) %>%
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      if(LoF_Homo != 0){
        Homo_name <- temp %>% 
          filter(IID == iid, Func == "LoF", Homozygous_M == T) %>% .$ID
        LoF_Homo_ID <- str_c(Homo_name, collapse = ";");Homo_name <- NULL
      } else {LoF_Homo_ID <- " "}
    }
    
    # CADD
    {
      for(value in c(0.03, 0.01)){
        CADD_Hetero <- temp %>% 
          filter(AF_eas <= value, PHRED >= CADD_score,IID == iid, Heterozygous_M == T) %>%
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
        
        if(CADD_Hetero != 0){
          Hetero_name <- temp %>% 
            filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Heterozygous_M == T) %>% .$ID
          CADD_Hetro_name <- str_c(Hetero_name, collapse = ";");Hetero_name <- NULL
        } else {CADD_Hetro_name <- " "}
        
        CADD_Homo <- temp %>% 
          filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Homozygous_M == T) %>%
          select(Clinical_Significance) %>%
          group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
        
        if(CADD_Homo != 0){
          Homo_name <- temp %>% 
            filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Homozygous_M == T) %>% .$ID
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
           CADD_MAF001_ALL = (CADD_MAF001_Hetero + (CADD_MAF001_Homo * 2))
           # ,Pathogenic_Hetro, Pathogenic_Hetro_ID,
           # Pathogenic_Homo, Pathogenic_Homo_ID,
           # LoF_Hetero, LoF_Hetro_ID, LoF_Homo, LoF_Homo_ID,
           # CADD_MAF003_Hetero, CADD_MAF003_Hetero_ID,
           # CADD_MAF003_Homo, CADD_MAF003_Homo_ID,
           # CADD_MAF001_Hetero, CADD_MAF001_Hetero_ID,
           # CADD_MAF001_Homo, CADD_MAF001_Homo_ID
    ) %>% return()
  }, mc.cores = detectCores() -1)
  
  return(sample_result)
} 
snp_table_per_sample_gene <- function(data_list, geneset_name, type, CADD_score = 20){
  index <- which(names(data_list) == geneset_name)
  temp <- data_list[[index]]
  temp$Clinical_Significance <- str_replace(temp$Clinical_Significance, ## missing resolve 0717
                                            pattern = "(^Pathogenic$)|^(Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic,_risk_factor$)",
                                            replacement = "Pathogenic_ALL")
  temp$Clinical_Significance[is.na(temp$Clinical_Significance)] <- "not_available"
  sample <- temp$IID %>% unique() %>% sort()
  
  # colnames set 
  name_change <- function(ID_DF, type){
    
    type <- paste0(type, "_Variant")
    col_n <- c()
    
    if(ncol(ID_DF) >= 1){
      for(index in 1:ncol(ID_DF)) {
        col_n <- c(col_n, paste0(type, index))
      }
      colnames(ID_DF) <- col_n  
      
      return(ID_DF)  
    } else return(ID_DF)
    
  }
  
  sample_result <- mclapply(X = sample, FUN = function(iid){
    # result <- list()
    PHENOTYPE <- temp %>% filter(IID == iid) %>% select(PHENOTYPE) %>% pull(1) %>% 
      unique() 
    PHENOTYPE <- ifelse(PHENOTYPE == 1, "control", "case")
    
    # Pathogenic ====
    Pathogenic <- temp %>% 
      filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
             IID == iid, Heterozygous != T) %>% 
      select(Clinical_Significance) %>%
      group_by_all() %>% count() %>% pull(2) %>% ifelse(length(.) == 0, 0, .)
    
    if(Pathogenic != 0){
      Pathogenic_name_He <- temp %>% 
        filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Heterozygous_M == T) %>% .$ID2 %>% as_tibble()
      Pathogenic_ID_He <- temp %>% 
        filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Heterozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
        bind_cols(., Pathogenic_name_He) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
        t() %>% as_tibble() %>% 
        name_change(ID_DF = ., type = "Pathogenic_Hetero")
      
      Pathogenic_name_Ho <- temp %>% 
        filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
               IID == iid, Homozygous_M == T) %>% .$ID2 %>% as_tibble()
      
      if(nrow(Pathogenic_name_Ho) == 0){
        Pathogenic_ID_Ho <- NULL
      } else{
        Pathogenic_ID_Ho <- temp %>% 
          filter(AF_eas <= 0.03, Clinical_Significance == "Pathogenic_ALL", 
                 IID == iid, Homozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
          bind_cols(., Pathogenic_name_Ho) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
          t() %>% as_tibble() %>% 
          name_change(ID_DF = ., type = "Pathogenic_Homo")
        # Homo mutation duplication!
        # Pathogenic_ID_Ho <- bind_cols(Pathogenic_ID_Ho, Pathogenic_ID_Ho)
      }
      
      Pathogenic_ID <- bind_cols(Pathogenic_ID_He, Pathogenic_ID_Ho)
      # Pathogenic_ID <- bind_cols(Pathogenic_ID_He, Pathogenic_ID_Ho) %>% 
      #   name_change(ID_DF = ., type = "Pathogenic")
      
    } else {Pathogenic_ID <- tibble(.rows = 1)}
    
    
    # LoF ====
    LoF <- temp %>% 
      filter(IID == iid, Func == "LoF", Heterozygous != T) %>%
      select(Clinical_Significance) %>%
      group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
    
    if(LoF!= 0){
      # Hetero mutation
      LoF_name_He <- temp %>% 
        filter(IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$ID2 %>% as_tibble()
      
      LoF_ID_He <- temp %>% 
        filter(IID == iid, Func == "LoF", Heterozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
        bind_cols(., LoF_name_He) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
        t() %>% as_tibble() %>% 
        name_change(ID_DF = ., type = "LoF_Hetero")
      
      # Homo mutation
      LoF_name_Ho <- temp %>% 
        filter(IID == iid, Func == "LoF", Homozygous_M == T) %>% .$ID2 %>% as_tibble()
      
      if(nrow(LoF_name_Ho) == 0){
        LoF_ID_Ho <- NULL
      } else {
        LoF_ID_Ho <- temp %>% 
          filter(IID == iid, Func == "LoF", Homozygous_M == T) %>% .$Gene.knownGene %>% as_tibble() %>% 
          bind_cols(., LoF_name_Ho) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
          t() %>% as_tibble() %>% 
          name_change(ID_DF = ., type = "LoF_Homo")
        # # Homo mutation duplication
        # LoF_ID_Ho <- bind_cols(LoF_ID_Ho, LoF_ID_Ho)
      }
      
      LoF_ID <- bind_cols(LoF_ID_He, LoF_ID_Ho)
      # LoF_ID <- bind_cols(LoF_ID_He, LoF_ID_Ho) %>% 
      #   name_change(ID_DF = ., type = "LoF")
      
      
    } else {
      LoF_ID <- tibble(.rows = 1)
    }
    
    # CADD ====
    for(value in c(0.03, 0.01)){
      CADD <- temp %>% 
        filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Heterozygous != T) %>% select(Clinical_Significance) %>%        
        group_by_all() %>% count() %>% pull(2) %>% sum() %>% ifelse(length(.) == 0, 0, .)
      
      if(CADD != 0){
        # Hetero Mutation Check!
        CADD_name_He <- temp %>% 
          filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Heterozygous_M == T) %>% 
          .$ID2 %>% as_tibble()
        CADD_ID_He <- temp %>% 
          filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Heterozygous_M == T) %>% 
          .$Gene.knownGene %>% as_tibble() %>% 
          bind_cols(., CADD_name_He) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
          t() %>% as_tibble() %>% 
          name_change(ID_DF = ., type = paste0("CADD_MAF", value, "_Hetero"))
        
        # Homo Mutation Check!
        CADD_name_Ho <- temp %>% 
          filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Homozygous_M == T) %>% 
          .$ID2 %>% as_tibble()
        
        if(nrow(CADD_name_Ho) == 0){
          CADD_ID_Ho <- NULL
          
        } else{
          CADD_ID_Ho <- temp %>% 
            filter(AF_eas <= value, PHRED >= CADD_score, IID == iid, Homozygous_M == T) %>% 
            .$Gene.knownGene %>% as_tibble() %>% 
            bind_cols(., CADD_name_Ho) %>% transmute(gene_rs = paste0(value...1, "(", value...2,")")) %>%
            t() %>% as_tibble() %>% 
            name_change(ID_DF = ., type = paste0("CADD_MAF", value, "_HOMO"))
          
          # Homo mutation duplication!! 
          # CADD_ID_Ho <- bind_cols(CADD_ID_Ho,CADD_ID_Ho)
        }
        
        CADD_ID <- bind_cols(CADD_ID_Ho, CADD_ID_He) 
        
      } else {
        CADD_ID <- tibble(.rows = 1)
        CADD_ID_Ho_verfy <- NULL
      }
      
      
      if(value == 0.03){
        CADD_MAF003_ID <- CADD_ID
        # CADD_MAF003_ID <- CADD_ID %>%
        #   name_change(ID_DF = ., type = "CADD_MAF003")
        
      } else {
        CADD_MAF001_ID <- CADD_ID 
        # CADD_MAF001_ID <- CADD_ID %>%
        #   name_change(ID_DF = ., type = "CADD_MAF001")
      }
      
    }
    ####### ====
    snp_table <- switch(type, "Pathogenic" = Pathogenic_ID,
                        "LoF" = bind_cols(Pathogenic_ID, LoF_ID),
                        "CADD_MAF003" = bind_cols(Pathogenic_ID,CADD_MAF003_ID),
                        "CADD_MAF001" = bind_cols(Pathogenic_ID, CADD_MAF001_ID),
                        "ALL" = bind_cols(Pathogenic_ID, CADD_MAF003_ID, CADD_MAF001_ID, LoF_ID) %>% as_tibble() # only CADD 003
    )
    
    tibble(Subject_ID = iid, PHENOTYPE) %>% bind_cols(snp_table) %>% return()
  }, mc.cores = detectCores() -1)
  
  return(sample_result) 
}
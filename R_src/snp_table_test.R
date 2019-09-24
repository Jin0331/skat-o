# source & clinvar load ====
rm(list=ls());gc()
source("/home/jinoo/skat-o/R_src/SKAT_SNPTABLE_F_ver_04.R")

library_load()
# clinvar <- clinvar_load()
clinvar <- clinvar_db()
clinvar <- clinvar_load()

geneset <- geneset_load_SKAT()


# TABLE 2 ====
  geneset <- geneset_load_SKAT()
  # 3 = O2, 6 = PARK, 8 = NONPARK
  gene <- 6
  geneset_ <- geneset[[1]][[gene]][!is.na(geneset[[1]][[gene]])]
  
  # row_all_variants
  
  for(data_name in c("WES_merge", "PPMI", "NeuroX")){
    variant <- list()
    fix <- fix_load(data_name = data_name, type = "ROW")
    
    for(index in 1:length(geneset_)){
      variant[[index]] <- subset(fix, subset = ((Gene.knownGene == geneset_[index])))
    }
    bind_rows(variant) %>% nrow() %>% print()
  }
  
  # WES_merge
  table2_WES_merge <- table2_calc(QC_fix = fix_load("WES_merge", type = "QC"), geneset = geneset_)
  # PPMI
  table2_PPMI <- table2_calc(QC_fix = fix_load("PPMI", type = "QC"), geneset = geneset_)
  # NeuroX
  table2_NeuroX <- table2_calc(QC_fix = fix_load("NeuroX", type = "QC"), geneset = geneset_)
  
# TABLE 3 WES all variant, CADD score upload & #CHROM ====

  WES_table <- snp_table_ALL(data_name = "WES_merge", index_ = c(3, 6, 8), clinvar = clinvar)
  PPMI_table <- snp_table_ALL(data_name = "PPMI", index_ = c(3, 6, 8), clinvar = clinvar)
  # NeuroX_table <- snp_table_ALL(index = c(3), data_name = "NeuroX", clinvar = clinvar)
  
  write_delim(x = bind_rows(WES_table), path = "WES_table_0830.txt", delim = "\t")
  write_delim(x = bind_rows(PPMI_table), path = "PPMI_table_0830.txt", delim = "\t")


# TABLE 5 6 preprocessing per individual ====
  wes_result <- snp_table(data_name = "WES_merge", index_ = c(3,6,8), clinvar = clinvar) 
  ppmi_result <- snp_table(data_name = "PPMI", index_ = c(3,6,8), clinvar = clinvar) 

# TABLE 5 ====
  
  for(geneset_name in c("O2", "O2_PARK", "O2_NONPARK")){
    table3_ver_1(data_list = wes_result, type = "Pathogenic", geneset_name = geneset_name) %>% 
      bind_cols(., table3_ver_1(data_list = wes_result, type = "LoF", geneset_name = geneset_name)) %>%
      bind_cols(., table3_ver_1(data_list = wes_result, type = "CADD_MAF", geneset_name = geneset_name, MAF_value = 0.03)) %>%
      bind_cols(., table3_ver_1(data_list = wes_result, type = "CADD_MAF", geneset_name = geneset_name, MAF_value = 0.01)) %>%
      bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname") %>% 
      write_delim(x = .,path = paste0("table5_result_WES_merge_",geneset_name,".txt"), delim = "\t")
  }
  
  for(geneset_name in c("O2", "O2_PARK", "O2_NONPARK")){
    table3_ver_1(data_list = ppmi_result, type = "Pathogenic", geneset_name = geneset_name, data_name = "PPMI") %>% 
      bind_cols(., table3_ver_1(data_list = ppmi_result, type = "LoF", geneset_name = geneset_name, data_name = "PPMI")) %>%
      bind_cols(., table3_ver_1(data_list = ppmi_result, type = "CADD_MAF", geneset_name = geneset_name, data_name = "PPMI", MAF_value = 0.03)) %>%
      bind_cols(., table3_ver_1(data_list = ppmi_result, type = "CADD_MAF", geneset_name = geneset_name, data_name = "PPMI", MAF_value = 0.01)) %>%
      bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname") %>% 
      write_delim(x = .,path = paste0("table5_result_PPMI_",geneset_name,".txt"), delim = "\t")
  }

# TABLE 6 ====
  # sample deleteriou
  data_name <- "WES_merge"
  FAM <- SKAT::Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/set2/WES_merge_0821.fam"),
                                  File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/set2/WES_merge_0821.cov"), Is.binary = FALSE) %>% as_tibble() %>%
    select(IID, AGE) %>% rename(Subject_ID = IID) %>% mutate(Subject_ID = as.character(Subject_ID))
  
  for(geneset_name in c("O2","O2_PARK","O2_NONPARK")){
    table3_ver_3_sample(data_list = wes_result, geneset_name = geneset_name) %>% bind_rows() %>%
      left_join(x = ., y = FAM, by = "Subject_ID") %>% 
      write_delim(x =., path = paste0(geneset_name,"_sample_WES_snv_.txt"), delim = "\t")
  }
  
  ## sample_gene
  O2_sample <- table3_ver_3_sample_gene(data_list = wes_result, geneset_name = "O2", type = "ALL") %>%
    bind_rows()
  O2_sample[is.na(O2_sample)] <- " "
  O2_sample_result <- O2_sample[, 3:ncol(O2_sample)] %>%
      select(sort(current_vars())) %>% bind_cols(O2_sample[, 1:2],.)
  write_delim(x = O2_sample_result, path = "WES_sample_O2_gene_0829.txt", delim = "\t")
  
  
  #### PPMI
  data_name <- "PPMI"
  FAM <- SKAT::Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/set2/PPMI_0821.fam"),
                                  File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/set2/PPMI_0821.cov"), Is.binary = FALSE) %>% as_tibble() %>%
    select(IID, AGE) %>% rename(Subject_ID = IID) %>% mutate(Subject_ID = as.character(Subject_ID))
  
  for(geneset_name in c("O2","O2_PARK","O2_NONPARK")){
    table3_ver_3_sample(data_list = ppmi_result, geneset_name = geneset_name) %>% bind_rows() %>%
      left_join(x = ., y = FAM, by = "Subject_ID") %>% 
      write_delim(x =., path = paste0(geneset_name,"_sample_PPMI_snv_.txt"), delim = "\t")
  }
  
  ## sample_gene
  O2_sample <- table3_ver_3_sample_gene(data_list = ppmi_result, geneset_name = "O2", type = "ALL") %>%
    bind_rows()
  O2_sample[is.na(O2_sample)] <- " "
  O2_sample_result <- O2_sample[, 3:ncol(O2_sample)] %>%
    select(sort(current_vars())) %>% bind_cols(O2_sample[, 1:2],.)
  write_delim(x = O2_sample_result, path = "PPMI_sample_O2_gene_0829.txt", delim = "\t")
  
  
  
# NOT RUN ====
  ## gene
  {
    O2_PARK_sup <- table3_ver_2_gene(WES_table = WES_table, geneset_name = "O2_PARK") %>% bind_rows()
    O2_NONPARK_sup <- table3_ver_2_gene(WES_table = WES_table, geneset_name = "O2_NONPARK") %>% bind_rows()
    
    write_delim(x = bind_rows(O2_PARK_sup, O2_NONPARK_sup), path = "Suppl_summary_of_deleter_snvs.txt", delim = "\t")
  }







## SNV SNV
{
  snv_snv <- bind_rows(WES_table) %>% select(CHROM, POS, ID, REF, ALT, Gene.knownGene, 
                                                AAChange.knownGene, MAF, Clinical_Significance, geneset, Func, CADD13_PHRED)
    # mutate(PHENOTYPE = ifelse(PHENOTYPE == 1, "control", "case"))
  
  write_delim(x = snv_snv, path = "Suppl_snv_snv.txt", delim = "\t")
}
{
  # ipdgc_table3_003 <- map(.x = ipdgc_result, .f = function(x){
  #   temp <- as.tibble(x) %>% arrange(IID) %>% 
  #     filter(MAF < 0.03) %>% 
  #     str_replace(.$Clinical_Significance, pattern = "(^Pathogenic$)|(^Likely_pathogenic$)|(^Pathogenic/Likely_pathogenic$)",
  #                 replacement = "Pathogenic_ALL") %>%
  #     str_replace_na(.$Clinical_Significance, replacement = "not_available")
  #     table3_CLSIG() %>% bind_rows()
  #   temp[is.na(temp)] <- 0
  #   
  #   control <- temp %>% filter(PHENOTYPE == 1) %>% table3_calc()
  #   case <- temp %>% filter(PHENOTYPE == 2) %>% table3_calc()
  #   
  #   return(list(control = control, case = case))
  # }) %T>% table3_write(table = ., data_name = "IPDGC", 0.03)
  # ipdgc_table3_001 <- map(.x = ipdgc_result, .f = function(x){
  #   temp <- x %>% filter(AnnotationClinvar == "TRUE") %>% mutate(Clinical_Significance = as.factor(Clinical_Significance)) %>%
  #     arrange(IID) %>% filter(MAF < 0.01) %>% select(IID, Clinical_Significance, PHENOTYPE) %>% as_tibble() %>%
  #     table3_CLSIG() %>% bind_rows()
  #   temp[is.na(temp)] <- 0
  #   
  #   control <- temp %>% filter(PHENOTYPE == 1) %>% table3_calc()
  #   case <- temp %>% filter(PHENOTYPE == 2) %>% table3_calc()
  #   
  #   return(list(control = control, case = case))
  # }) %T>% table3_write(table = ., data_name = "IPDGC", 0.01)
  
  {
    # neurox_table3_003 <- map(.x = neuroX_result, .f = function(x){
    #   temp <- x %>% filter(AnnotationClinvar == "TRUE") %>% mutate(Clinical_Significance = as.factor(Clinical_Significance)) %>%
    #     arrange(IID) %>% filter(MAF < 0.03) %>% select(IID, Clinical_Significance, PHENOTYPE) %>% as_tibble() %>%
    #     table3_CLSIG() %>% bind_rows()
    #   temp[is.na(temp)] <- 0
    #   
    #   control <- temp %>% filter(PHENOTYPE == 1) %>% table3_calc()
    #   case <- temp %>% filter(PHENOTYPE == 2) %>% table3_calc()
    #   
    #   return(list(control = control, case = case))
    # }) %T>% table3_write(table = ., data_name = "NeuroX", 0.03)
    # neurox_table3_001 <- map(.x = neuroX_result, .f = function(x){
    #   temp <- x %>% filter(AnnotationClinvar == "TRUE") %>% mutate(Clinical_Significance = as.factor(Clinical_Significance)) %>%
    #     arrange(IID) %>% filter(MAF < 0.01) %>% select(IID, Clinical_Significance, PHENOTYPE) %>% as_tibble() %>%
    #     table3_CLSIG() %>% bind_rows()
    #   temp[is.na(temp)] <- 0
    #   
    #   control <- temp %>% filter(PHENOTYPE == 1) %>% table3_calc()
    #   case <- temp %>% filter(PHENOTYPE == 2) %>% table3_calc()
    #   
    #   return(list(control = control, case = case))
    # }) %T>% table3_write(table = ., data_name = "NeuroX", 0.01)
  }
}

rm(list=ls());gc()
# source("/home/jinoo/skat-o/snp_table.R")
source("/home/jinoo/skat-o/snp_table_0716.R")

library_load()
clinvar <- clinvar_load()

# ## TABLE 1 & 2
{
  
  
  
  # geneset load
  geneset <- geneset_load()
  
  # case control
  FAM <- data.table::fread(file = "/home/jinoo/skat-o/SKAT_data/IPDGC.fam", header = F)
  BIM <- data.table::fread(file = "/home/jinoo/skat-o/SKAT_data/IPDGC.bim", header = F)
  
  FAM <- data.table::fread(file = "/home/jinoo/skat-o/SKAT_data/NeuroX.fam", header = F)
  BIM <- data.table::fread(file = "/home/jinoo/skat-o/SKAT_data/NeuroX.bim", header = F)
  
  FAM %>% filter(V6 == 2) %>% nrow()
  
  target_gene <- c()
  
  for(gene in index)
    target_gene <- c(target_gene, geneset[[1]][[gene]][!is.na(geneset[[1]][[gene]])])
  variant_gene <- unique.default(target_gene)
  
  IPDGC_gene <- NULL;NeuroX_gene <- NULL;temp <- NULL
  for(i in 1:length(variant_gene)){
    temp <- subset.data.frame(IPDGC_FIX,
                              subset = (Gene.knownGene %in% variant_gene[i]))
    IPDGC_gene <- bind_rows(IPDGC_gene, temp)
  }
  for(i in 1:length(variant_gene)){
    temp <- subset.data.frame(NeuroX_Fix,
                              subset = (Gene.knownGene %in% variant_gene[i]))
    NeuroX_gene <- bind_rows(NeuroX_gene, temp)
  }
  
  
  
  IPDGC_gene %>% filter(ExonicFunc.knownGene == "nonsynonymous_SNV"
                        | ExonicFunc.knownGene == "stopgain"
                        | ExonicFunc.knownGene ==  "stoploss"
                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                        | ExonicFunc.knownGene ==  "splicing") %>% nrow()
  
  IPDGC_gene %>% filter(ExonicFunc.knownGene == "stopgain"
                        | ExonicFunc.knownGene ==  "stoploss"
                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                        | ExonicFunc.knownGene ==  "splicing")  %>%
    filter(MAF < 0.03) %>%
    nrow()
  
}

## TABLE 2
{
  test_fix <- WES_merge_fix <- fix_load("IPDGC", type = "ROW") # row fix ㅅㅜㅈㅓㅇ 
  test_fix <- WES_QC_fix <- fix_load("IPDGC", type = "QC")
  test_fix <- NeuroX_QC_fix <- fix_load("NeuroX", type = "ROW")
  test_fix <- NeuroX_QC_fix <- fix_load("NeuroX", type = "QC")
  geneset <- geneset_load()
}

# 3 = O2, 6 = PARK, 8 = NONPARK, 22 = CHD
{
  gene <- 3
  geneset_ <- geneset[[1]][[gene]][!is.na(geneset[[1]][[gene]])]
  variant <- list()
  
  
  # row_all_variants
  for(i in 1:length(geneset_)){
    variant[[i]] <- subset(test_fix, subset = ((Gene.knownGene == geneset_[i])))
  }
  
  bind_rows(variant) %>% nrow()
  
  ## Parkinson geneset variant
  for(i in 1:length(geneset_)){
    variant[[i]] <- subset(test_fix, subset = ((Gene.knownGene == geneset_[i] & Func.knownGene == "exonic")),
                           select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED","MAF"))
  }
  variant <- bind_rows(variant)
  
  ### 1. nonsynonymous geneset
  print(paste(gene, "non-syn", sep = " "))
  nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"
                                             | ExonicFunc.knownGene == "stopgain"
                                             | ExonicFunc.knownGene ==  "stoploss"
                                             | ExonicFunc.knownGene ==  "frameshift_deletion"
                                             | ExonicFunc.knownGene ==  "frameshift_insertion"
                                             | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                             | ExonicFunc.knownGene ==  "splicing") ) %>%
    nrow() %>% print()
  
  
  # ### 2. CADD > 12.37 variant
  # print(paste(gene, "cadd", sep = " "))
  # cadd <- subset(variant, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
  #                                    | ExonicFunc.knownGene == "stopgain"
  #                                    | ExonicFunc.knownGene ==  "stoploss"
  #                                    | ExonicFunc.knownGene ==  "frameshift_deletion"
  #                                    | ExonicFunc.knownGene ==  "frameshift_insertion"
  #                                    | ExonicFunc.knownGene ==  "frameshift_block_substitution"
  #                                    | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37) %>%
  #   nrow() %>% print()
  
  ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
  print(paste(gene, "lof", sep = " "))
  
  lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                    | ExonicFunc.knownGene ==  "stoploss"
                                    | ExonicFunc.knownGene ==  "frameshift_deletion"
                                    | ExonicFunc.knownGene ==  "frameshift_insertion"
                                    | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                    | ExonicFunc.knownGene ==  "splicing") ) %>%
    nrow() %>% print()
}

# TABLE 3 WES all variant, CADD score upload & #CHROM
{
  ipdgc_variant <- snp_table_WES(index_ = c(3,6,8,22), clinvar = clinvar) %>% bind_rows()
  neuroX_variant <- snp_table_WES(index_ = c(3,6,8,22), clinvar = clinvar, data_name = "NeuroX") %>% bind_rows()
  fwrite(x = select(ipdgc_variant, -CADD13_PHRED), file = "ipdgc_variant.vcf",sep = "\t", col.names = T, row.names = F)
  fwrite(x = select(neuroX_variant, -CADD13_PHRED), file = "neuroX_variant.vcf",sep = "\t", col.names = T, row.names = F)
  system("gzip ipdgc_variant.vcf")
  system("gzip neuroX_variant.vcf")
  
  CADD_annotation <- fread(file = "NeuroX_cadd.tsv.gz") %>%
    rename(CHROM = `#CHROM`, POS = POS, REF = REF, ALT = ALT, CADD13 = RawScore, CADD13_PHRED = PHRED) %>% 
    # mutate(CHROM = as.character(CHROM)) %>%
    left_join(x = neuroX_variant, y = ., by = c("CHROM","POS","REF","ALT")) %T>%
    fwrite(x = ., file = "NeuroX_all_variant.txt", col.names = T, row.names = F, sep = "\t")
}

# TABLE 3 preprocessing per individual
WES_table <- snp_table_ALL(index_ = c(3,6,8), clinvar = clinvar)
# NeuroX_table <- snp_table_ALL(index_ = c(3,6,8), clinvar = clinvar, data_name = "NeuroX")
ipdgc_result <- snp_table(data_name = "IPDGC", index_ = c(3,6,8), clinvar = clinvar)

# TABLE 3
## Pathogenic, LoF, CADD_MAF\
{
  table3_result_O2 <- table3_ver_1(data_list = ipdgc_result, type = "Pathogenic", geneset_name = "O2") %>% 
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "LoF", geneset_name = "O2")) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2", MAF_value = 0.03)) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2", MAF_value = 0.01)) %>%
    bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname")
  
  table3_result_O2_PARK <- table3_ver_1(data_list = ipdgc_result, type = "Pathogenic", geneset_name = "O2_PARK") %>% 
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "LoF", geneset_name = "O2_PARK")) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2_PARK", MAF_value = 0.03)) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2_PARK", MAF_value = 0.01)) %>%
    bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname")
  
  table3_result_O2_NONPARK <- table3_ver_1(data_list = ipdgc_result, type = "Pathogenic", geneset_name = "O2_NONPARK") %>% 
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "LoF", geneset_name = "O2_NONPARK")) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2_NONPARK", MAF_value = 0.03)) %>%
    bind_cols(., table3_ver_1(data_list = ipdgc_result, type = "CADD_MAF", geneset_name = "O2_NONPARK", MAF_value = 0.01)) %>%
    bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname")
  
  write_delim(x = table3_result_O2, path = "table3_result_O2.txt", delim = "\t")
  write_delim(x = table3_result_O2_NONPARK, path = "table3_result_O2_NONPARK.txt", delim = "\t")
  write_delim(x = table3_result_O2_PARK, path = "table3_result_O2_PARK.txt", delim = "\t")
}

# Suppl 
## gene
{
  O2_PARK_sup <- table3_ver_2_gene(WES_table = WES_table, geneset_name = "O2_PARK") %>% bind_rows()
  O2_NONPARK_sup <- table3_ver_2_gene(WES_table = WES_table, geneset_name = "O2_NONPARK") %>% bind_rows()
  
  write_delim(x = bind_rows(O2_PARK_sup, O2_NONPARK_sup), path = "Suppl_summary_of_deleter_snvs.txt", delim = "\t")
}

## sample
{
  data_name <- "IPDGC"
  FAM <- SKAT::Read_Plink_FAM_Cov(Filename = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".fam"),
                            File_Cov = paste0("/home/jinoo/skat-o/SKAT_data/",data_name,".cov"), Is.binary = FALSE) %>% as_tibble() %>%
    select(IID, AGE) %>% rename(Subject_ID = IID) %>% mutate(Subject_ID = as.character(Subject_ID))
  
  O2_sample <- table3_ver_3_sample(data_list = ipdgc_result, geneset_name = "O2") %>% bind_rows() %>%
    left_join(x = ., y = FAM, by = "Subject_ID")
  O2_PARK_sample <- table3_ver_3_sample(data_list = ipdgc_result, geneset_name = "O2_PARK") %>% bind_rows() %>%
    left_join(x = ., y = FAM, by = "Subject_ID")
  O2_NONPARK_sample <- table3_ver_3_sample(data_list = ipdgc_result, geneset_name = "O2_NONPARK") %>% bind_rows() %>%
    left_join(x = ., y = FAM, by = "Subject_ID")
  
  write_delim(x =O2_sample, path = "O2_sample_snv_.txt", delim = "\t")
  write_delim(x =O2_PARK_sample, path = "O2_PARK_sample_snv_.txt", delim = "\t")
  write_delim(x =O2_NONPARK_sample, path = "O2_NONPARK_snv_.txt", delim = "\t")
}

## sample_gene
{
  O2_sample <- table3_ver_3_sample_gene(data_list = ipdgc_result, geneset_name = "O2", type = "CADD_MAF003") %>%
    bind_rows()
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

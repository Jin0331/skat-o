# source & clinvar load ====
source("/home/dblab/skat-o/test_FPD/snp_table.R")
library_load()
plink_path <- "/home/dblab/QC_dir/QC/plink/FPD_missing_skat_QC"
clinvar <- clinvar_load(clinvar_path = "/home/dblab/temp/temp/clinvar_20210110.vcf")
geneset <- geneset_KFPD()
fix <- fix_load_KFPD()

#  All variant, CADD score upload & #CHROM ====
KFPD_table <- snp_table_variant(data_name = "KoreanFPD", index_ = 1:3, fix = fix, geneset = geneset, 
                                clinvar = clinvar, plink_path = plink_path)
write_delim(x = bind_rows(KFPD_table), path = "KoreanFPD_all_variant_info.txt", delim = "\t", na = "-")


# snp_table sample preprocessing per individual ====
KFPD_data_list <- snp_table_preprocessing(data_name = "KoreanFPD", index_ = 1:3, clinvar = clinvar, fix = fix, geneset = geneset,
                            plink_path = plink_path)

# snp_table_sample_count ====
for(geneset_name in names(KFPD_data_list)){
  snp_table_sample_count(data_list = KFPD_data_list, type = "Pathogenic", geneset_name = geneset_name) %>% 
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "LoF", geneset_name = geneset_name)) %>%
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "CADD_MAF", geneset_name = geneset_name, MAF_value = 0.03)) %>%
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "CADD_MAF", geneset_name = geneset_name, MAF_value = 0.01)) %>%
    bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname") %>% 
    write_delim(x = .,path = paste0("KoreanKFPD_sample_count_",geneset_name,".txt"), delim = "\t")
}

# sample deleteriou ====
# sample count
data_name <- "KoreanFPD"
FAM <- SKAT::Read_Plink_FAM_Cov(Filename = paste0(plink_path, ".fam"),
                                File_Cov = paste0(plink_path, ".cov"), Is.binary = FALSE) %>% as_tibble() %>%
  select(IID, AGE) %>% rename(Subject_ID = IID) %>% mutate(Subject_ID = as.character(Subject_ID))

temp_list <- list()
for(geneset_name in names(KFPD_data_list)){
  temp_list[[`geneset_name`]] <- snp_table_per_sample(data_list = KFPD_data_list, geneset_name = geneset_name) %>% bind_rows() %>%
    left_join(x = ., y = FAM, by = "Subject_ID") %>% 
    rename_at(.vars = 4:7, .funs = ~ paste(., geneset_name, sep = "_"))
    # write_delim(x =., path = paste0(geneset_name,"_sample_KFPD_snv_.txt"), delim = "\t")
}
temp_list <- temp_list[[1]] %>% inner_join(x = ., y = temp_list[[2]], by = c("Subject_ID", "SEX", "AGE", "PHENOTYPE")) %>% 
  inner_join(x = ., y = temp_list[[3]], by = c("Subject_ID", "SEX", "AGE", "PHENOTYPE"))


## sample gene
KFPD_sample <- snp_table_per_sample_gene(data_list = KFPD_data_list, geneset_name = names(KFPD_data_list)[3], type = "ALL") %>%
  bind_rows()
KFPD_sample[is.na(KFPD_sample)] <- " "
KFPD_sample_result <- KFPD_sample[, 3:ncol(KFPD_sample)] %>%
  select(sort(current_vars())) %>% bind_cols(KFPD_sample[, 1:2],.)

KPD_sample_gene <- left_join(x = temp_list, y = KFPD_sample_result, by = c("Subject_ID", "PHENOTYPE"))
KPD_sample_gene[is.na(KPD_sample_gene)] <- " "
write_delim(x = KPD_sample_gene, path = "KFPD_sample_gene.txt", delim = "\t")

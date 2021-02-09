# source & clinvar load ====
source("case_only_snptable_function.R")
library_load()
plink_path <- "/home/dblab/temp/case_only_snptable/data/all_case_105"
plink_path <- "/home/dblab/temp/case_only_snptable/data/exclude_case"
file_name <- "all_case_105"
# plink_path <- "/home/dblab/QC_dir/exclude_case/set1/exclude_case_id_phenosex"

clinvar <- clinvar_load(clinvar_path = "/home/dblab/temp/temp/clinvar_20210110.vcf")
geneset <- geneset_KFPD(path = "https://s3.us-west-2.amazonaws.com/secure.notion-static.com/7b12f84e-e41f-4956-babf-4df080ea1aae/input_gene.txt?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAT73L2G45O3KS52Y5%2F20210208%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20210208T073632Z&X-Amz-Expires=86400&X-Amz-Signature=87fc9bcca32e96979f876e2e13070e29993656e43e88ef15dedec110b4a27290&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22input_gene.txt%22")
fix <- fix_load_KFPD(path = "/home/dblab/temp/case_only_snptable/data/all_case_106_fix_anno.txt",
                     plink_path = plink_path)

dir.create("test")
setwd("test")
#  All variant, CADD score upload & #CHROM ====
KFPD_table <- snp_table_variant(data_name = file_name, index_ = 1:3, fix = fix, geneset = geneset, 
                                clinvar = clinvar, plink_path = plink_path)
write_delim(x = bind_rows(KFPD_table), path = paste0("../", file_name,"_snp_table.txt"), delim = "\t", na = "-")


# snp_table sample preprocessing per individual ====
KFPD_data_list <- snp_table_preprocessing(data_name = file_name, index_ = 1:3, clinvar = clinvar, fix = fix, geneset = geneset,
                                          plink_path = plink_path)

# snp_table_sample_count ====
case_number <- 105
control_number <- 0
for(geneset_name in names(KFPD_data_list)){
  snp_table_sample_count(data_list = KFPD_data_list, type = "Pathogenic", 
                         geneset_name = geneset_name, case_number = case_number, control_number = control_number) %>% 
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "LoF", 
                                        geneset_name = geneset_name, case_number = case_number, control_number = control_number)) %>%
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "CADD_MAF", 
                                        geneset_name = geneset_name, MAF_value = 0.03, case_number = case_number, control_number = control_number)) %>%
    bind_cols(., snp_table_sample_count(data_list = KFPD_data_list, type = "CADD_MAF", 
                                        geneset_name = geneset_name, MAF_value = 0.01, case_number = case_number, control_number = control_number)) %>%
    bind_cols(tibble(rowname = c("0","1","2","3","4","5","6","7","8+")), .) %>% column_to_rownames(., var = "rowname") %>% 
    write_delim(x = .,path = paste0(file_name,"_sample_count_",geneset_name,".txt"), delim = "\t")
}

# sample deleteriou ====
# sample count

temp_list <- list()
for(geneset_name in names(KFPD_data_list)){
  temp_list[[`geneset_name`]] <- snp_table_per_sample(data_list = KFPD_data_list, geneset_name = geneset_name) %>% bind_rows() %>%
    rename_at(.vars = 4:7, .funs = ~ paste(., geneset_name, sep = "_"))
}

temp_list <- temp_list[[1]] %>% inner_join(x = ., y = temp_list[[2]], by = c("Subject_ID", "SEX", "PHENOTYPE")) %>% 
  inner_join(x = ., y = temp_list[[3]], by = c("Subject_ID", "SEX", "PHENOTYPE"))


## sample gene
KFPD_sample <- snp_table_per_sample_gene(data_list = KFPD_data_list, geneset_name = names(KFPD_data_list)[3], type = "ALL") %>%
  bind_rows()
KFPD_sample[is.na(KFPD_sample)] <- " "
KFPD_sample_result <- KFPD_sample[, 3:ncol(KFPD_sample)] %>%
  select(sort(current_vars())) %>% bind_cols(KFPD_sample[, 1:2],.)

KPD_sample_gene <- left_join(x = temp_list, y = KFPD_sample_result, by = c("Subject_ID", "PHENOTYPE"))
KPD_sample_gene[is.na(KPD_sample_gene)] <- " "
write_delim(x = KPD_sample_gene, path = paste0(file_name, "_sample_gene.txt"), delim = "\t")

source("/home/jinoo/skat-o/R_src/SKAT_SNPTABLE_F_ver_01.R")

geneset <- geneset_load_SKAT()
for(dataset in c("IPDGC","NeuroX")){
  fix <- fix_load_SKAT(dataset)
  
  for(cadd_score in c(20)){
  # geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(3,4,6,7,8,9),
                # CADD_score = cadd_score)
  gene_setid(geneset = geneset, fix, data_name = dataset, index = c(3), CADD_score = cadd_score)
  # run_skat_all_cov(data_name = dataset, add_name = as.character(cadd_score))
  run_skat_all_cov(data_name = dataset, flag = "per_gene",add_name = as.character(cadd_score))
  MAC_calculation(data_name = dataset, add_name = "0806")
  
  # run_skat_all_cov(data_name = dataset, flag = "gene")
  # run_skat_all_common_rare_cov(data_name = dataset)
  }
}


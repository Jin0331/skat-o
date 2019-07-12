source("/home/jinoo/skat-o/skat_function.R")

for(dataset in c("IPDGC")){
  fix <- fix_load(dataset)
  geneset <- geneset_load()
  geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(3),
                CADD_score = 12.37)
  # gene_setid(geneset = geneset, fix, data_name = dataset, index = c(4,9), CADD_score = 23)
  
  run_skat_all_cov(data_name = dataset)
  # run_skat_all_cov(data_name = dataset, flag = "gene")
  # run_skat_all_common_rare_cov(data_name = dataset)
}

source("/home/jinoo/skat-o/skat_function.R")

for(dataset in c("IPDGC","NeuroX")){
  fix <- fix_load(dataset)
  geneset <- geneset_load()
  # geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(2, 14:15, 17:19, 23))
  gene_setid(geneset = geneset, fix, data_name = dataset, index = 26)
  
  # run_skat_all_cov(data_name = dataset)
  run_skat_all_cov(data_name = dataset, flag = "gene")
  # run_skat_all_common_rare_cov(data_name = dataset)
}

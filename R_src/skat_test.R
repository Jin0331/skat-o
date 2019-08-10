# source & init var =====
source("/home/jinoo/skat-o/R_src/SKAT_SNPTABLE_F_ver_01.R")
library_load()
#c(3,4,6,7,8,9)
geneset <- geneset_load_SKAT()
cadd_score <- 20
add_name <- "cadd_20"

# SKAT-O test geneset & per gene =====
for(dataset in c("IPDGC","NeuroX")){
  if(dataset == "IPDGC"){
    for(set_IPDGC in c("set1","set2","set3")){
      fix <- fix_load_SKAT(dataset, set_IPDGC = set_IPDGC)
      geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(3,4,6,7,8,9),
                    CADD_score = cadd_score)
      run_skat_all_cov(data_name = dataset, add_name = as.character(cadd_score), set_IPDGC = set_IPDGC) %>% 
        bind_rows() %>% p_adjuste_cal() %>% 
        fwrite(x = ., file = paste0(dataset,"_result_cov_all_", add_name,"_", set_IPDGC,".txt"), 
               col.names = T, row.names = F, sep = "\t")
      
      file.remove(list.files()[str_detect(list.files(), pattern = "SetID")])}
    
  } else {
    fix <- fix_load_SKAT(dataset)
    geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(3,4,6,7,8,9),
                  CADD_score = cadd_score)
    run_skat_all_cov(data_name = dataset, add_name = as.character(cadd_score)) %>% 
      p_adjuste_cal() %>% 
      fwrite(x = ., file = paste0(dataset,"_result_cov_all_", add_name,".txt"), 
             col.names = T, row.names = F, sep = "\t")
    file.remove(list.files()[str_detect(list.files(), pattern = "SetID")])
  }
}

# NOT RUN ====
# gene_setid(geneset = geneset, fix, data_name = dataset, index = c(3), CADD_score = cadd_score)
# run_skat_all_cov(data_name = dataset, flag = "per_gene",add_name = as.character(cadd_score))
# MAC_calculation(data_name = dataset, add_name = "0806")
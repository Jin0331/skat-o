# source & init var =====
source("/home/jinoo/skat-o/R_src/SKAT_SNPTABLE_F_ver_04.R")
# source("/home/jinoo/skat-o/R_src/test_0917.R")
library_load()
#c(3,4,6,7,8,9)
#c(27,28,29,30)

#c(2,3,4,6,7,8,9,30) # final geneset 0830
geneset <- geneset_load_SKAT()


if(!dir.exists("test")){
  dir.create("test")
  setwd("test/")
}

setwd("test/")

# SKAT-O test geneset & per gene =====
for(dataset in c("WES_merge")){
  
  cadd_score <- 20
  resampling <- 0
  add_name <- "0919"
  
  for(set_WES in c("set2")){
    if(getwd() != paste0("/home/jinoo/skat-o/SKAT_data/", set_WES,"/test"))
      setwd(paste0("/home/jinoo/skat-o/SKAT_data/", set_WES,"/test"))
    fix <- fix_load_SKAT(dataset, set_WES  = set_WES)
    geneset_setid(geneset_merge = geneset[[1]], geneset[[2]], fix, data_name = dataset, index_ = c(2,3,4,6,7,8,9,30), CADD_score = cadd_score)
    # gene_setid(geneset = geneset, fix, data_name = dataset, index = c(3), CADD_score = cadd_score)
    
    result <- run_skat_all_cov(data_name = dataset, add_name = as.character(cadd_score), set_WES = set_WES,
                               re = resampling) %>%
      bind_rows() %>% #p_adjuste_cal() %>%
      fwrite(x = ., file = paste0(dataset,"_SKATO_", add_name,"_", set_WES,".txt"),
             col.names = T, row.names = F, sep = "\t")

    # run_skat_all_cov(data_name = dataset, flag = "gene",add_name = as.character(cadd_score), set_WES = set_WES) %>%
    #    bind_rows() %>% #p_adjuste_cal() %>%
    #    fwrite(x = ., file = paste0(dataset,"_SKATO_per_gene", add_name,"_", set_WES,".txt"),
    #           col.names = T, row.names = F, sep = "\t")
    
    file.remove(list.files()[str_detect(list.files(), pattern = "SetID|bim|fam|bed")])
  }
} 



 
 # data_name <- "IPDGC"
# data_name <- "NeuroX"
# 
# fix <- fix_load_SKAT(data_name = data_name, set_IPDGC = "set1")
# gene_setid(geneset = geneset, fix, data_name = data_name, index = c(3), CADD_score = cadd_score)
# MAC_calculation(data_name = data_name, add_name = "0806", set_IPDGC = "set1")


# NOT RUN ====
# gene_setid(geneset = geneset, fix, data_name = dataset, index = c(3), CADD_score = cadd_score)
# run_skat_all_cov(data_name = dataset, flag = "per_gene",add_name = as.character(cadd_score))
# MAC_calculation(data_name = dataset, add_name = "0806")
# source & init var =====
source("/home/jinoo/skat-o/R_src/SKAT_SNPTABLE_F_ver_04.R")
# source("/home/jinoo/skat-o/R_src/test_0917.R")
library_load()
#c(3,4,6,7,8,9)
#c(27,28,29,30)

#c(2,3,4,6,7,8,9,30) # final geneset 0830
#geneset <- geneset_load_SKAT()
load("/home/jinoo/skat-o/SKAT_data/geneset_fix_clinvar/geneset_0805.RData")

if(!dir.exists("/home/jinoo/skat-o/test")){
  dir.create("/home/jinoo/skat-o/test")
  setwd("/home/jinoo/skat-o/test")
} else{
  setwd("/home/jinoo/skat-o/test")
}

# Meta-analysis ====
for(dataset in c("WES_merge", "PPMI","NeuroX")){
  cadd_score <- 20
  resampling <- 0
  add_name <- "0608_test_"
  
  for(set_WES in c("set2")){
    fix <- fix_load_SKAT(dataset, set_WES  = set_WES)
  
    # geneset test!
    geneset_setid(geneset_merge = geneset[[1]], col_name = geneset[[2]], test_fix = fix, 
                  data_name = dataset, index_ = c(3,6,8), CADD_score = cadd_score, geneset_type = "ALL")
    
    obj <- make_obj(data_name = dataset, add_name = as.character(cadd_score), set_WES = set_WES, re = resampling)
    N.Sample <- read.table(paste0("/home/jinoo/skat-o/SKAT_data/", set_WES, "/",dataset, "_0821.fam")) %>% 
      .[, 6] %>% length()
    
    Generate_Meta_Files(obj = obj, paste0("/home/jinoo/skat-o/SKAT_data/", set_WES, "/",dataset, "_0821.bed"),
                          paste0("/home/jinoo/skat-o/SKAT_data/", set_WES, "/",dataset, "_0821.bim"),
                          paste0(dataset,"_skat.SetID"),
                          paste0(dataset,".MSSD"),
                          paste0(dataset,".MInfo"),
                          N.Sample)
    
  }
}

file.remove(list.files()[str_detect(list.files(), pattern = "SetID|bim|fam|bed")])

File.Mat.vec <- c()
File.Info.vec <- c()
for(dataset in c("WES_merge", "PPMI","NeuroX")){
  File.Mat <- sprintf("./%s.MSSD", dataset)
  File.Info <- sprintf("./%s.MInfo", dataset)
  
  File.Mat.vec <- c(File.Mat.vec, File.Mat)
  File.Info.vec <- c(File.Info.vec, File.Info)
}

# heterogenuos
Cohort.Info<-Open_MSSD_File_2Read(File.Mat.vec, File.Info.vec)
maf_003_cw <- MetaSKAT_MSSD_ALL(Cohort.Info, method = "optimal", Group_Idx = NULL,
                                MAF.cutoff = 0.03, combined.weight = FALSE)
write_delim(x = maf_003_cw, path = "maf_003_cw.txt", delim = "\t", col_names = T)

maf_001_cw <- MetaSKAT_MSSD_ALL(Cohort.Info, method = "optimal", Group_Idx = NULL, 
                                MAF.cutoff = 0.01, combined.weight = FALSE)
write_delim(x = maf_001_cw, path = "maf_001_cw.txt", delim = "\t", col_names = T)

#WES, PPMI ---> SAME, NeuroX
maf_003_idx1 <- MetaSKAT_MSSD_ALL(Cohort.Info, method = "optimal", Group_Idx = c(1,1,2), 
                                  is.separate = TRUE, combined.weight = FALSE, MAF.cutoff = 0.03)
write_delim(x = maf_003_idx1, path = "maf_003_idx.txt", delim = "\t", col_names = T)

maf_001_idx1 <- MetaSKAT_MSSD_ALL(Cohort.Info, method = "optimal", Group_Idx = c(1,1,2), 
                                  is.separate = TRUE, combined.weight = FALSE, MAF.cutoff = 0.01)
write_delim(x = maf_001_idx1, path = "maf_001_idx.txt", delim = "\t", col_names = T)


# SKAT-O test geneset & per gene =====
# WES_merge, PPMI, NeuroX
for(dataset in c("WES_merge", "PPMI", "NeuroX")){
  cadd_score <- 20
  resampling <- 0
  add_name <- "excluding_GBA_LRRK2"
  
  for(set_WES in c("set2")){
    # if(getwd() != paste0("/home/jinoo/skat-o/SKAT_data/", set_WES,"/test"))
    #   setwd(paste0("/home/jinoo/skat-o/SKAT_data/", set_WES,"/test"))
    fix <- fix_load_SKAT(dataset, set_WES  = set_WES)
    
    # geneset test!   # c(2,3,6,8,23,24)
    geneset_setid(geneset_merge = geneset[[1]], 
                  col_name = geneset[[2]], 
                  test_fix = fix, 
                  data_name = dataset, 
                  index_ = c(4), 
                  CADD_score = cadd_score,
                  geneset_type = "ALL")
    
    result <- run_skat_all_cov(data_name = dataset, add_name = as.character(cadd_score), set_WES = set_WES, re = resampling) %>%
      bind_rows() %>% p_adjuste_cal() %>%
      fwrite(x = ., file = paste0(dataset,"_SKATO_", add_name,"_", set_WES,".txt"),
               col.names = T, row.names = F, sep = "\t")
      
    # per gene test!
    # gene_setid(geneset = geneset, fix, data_name = dataset, index = c(24), CADD_score = cadd_score)
    # run_skat_all_cov(data_name = dataset, flag = "gene",add_name = as.character(cadd_score), set_WES = set_WES) %>%
    #    bind_rows() %>% p_adjuste_cal() %>%
    #    fwrite(x = ., file = paste0(dataset,"_SKATO_per_gene", add_name,"_", set_WES,".txt"),
    #           col.names = T, row.names = F, sep = "\t")
    
    # file.remove(list.files()[str_detect(list.files(), pattern = "SetID|bim|fam|bed")])
  }
} 



 
 # data_name <- "IPDGC"
# data_name <- "NeuroX"
# 
# fix <- fix_load_SKAT(data_name = data_name, set_IPDGC = "set1")
# gene_setid(geneset = geneset, fix, data_name = data_name, index = c(3), CADD_score = cadd_score)
# MAC_calculation(data_name = data_name, add_name = "0806", set_IPDGC = "set1")



# power calc
WES_power <- skat_powerCalc(data_name = "WES_merge", region_type = "avg")
PPMI_power <- skat_powerCalc(data_name = "PPMI", region_type = "avg")
NeuroX_power <- skat_powerCalc(data_name = "NeuroX", region_type = "avg", set_seed = 400) %>% as_tibble() %>% 
  bind_cols(tibble(seed = 100))

WES_power_geneset <- skat_powerCalc(data_name = "WES_merge", region_type = "sum")
PPMI_power_geneset <- skat_powerCalc(data_name = "PPMI", region_type = "sum")
NeuroX_power_geneset <- skat_powerCalc(data_name = "NeuroX", region_type = "sum")


se_value <- sample(1:10000, 2345)
temp <- mclapply(se_value, function(value){
  skat_powerCalc(data_name = "NeuroX", region_type = "avg", set_seed = value) %>% as_tibble() %>% 
    bind_cols(tibble(seed = value)) %>% return()
}, mc.cores = 6)


# NOT RUN ====
# gene_setid(geneset = geneset, fix, data_name = dataset, index = c(3), CADD_score = cadd_score)
# run_skat_all_cov(data_name = dataset, flag = "per_gene",add_name = as.character(cadd_score))
# MAC_calculation(data_name = dataset, add_name = "0806")
library(vcfR);library(tidyverse);library(doParallel)
library(data.table)
# sample
{
  temp <- read.vcfR(file = "C:/Users/JINOO/Documents/ExAC.r0.3.1.sites.vep.vcf", verbose = T) 
  temp_tidy <- temp@fix %>% as_tibble() %>% select(-ID, -QUAL, -FILTER, -REF, -ALT)
  temp_info <- temp_tidy$INFO[1] %>% str_split(string = ., pattern = ";", simplify = T) %>% str_split(string = ., pattern = "=", simplify = T) %>%
    as_tibble() %>% mutate(V2 = as.numeric(V2)) %>% subset.data.frame(., subset = (V1 == "AC_Adj" | V1 == "AN_Adj"))
  AC_AF <- tibble(AF_adj = temp_info$V2[1]/temp_info$V2[2])
}

# ExAC MAF AC_Adj / AN_Adj
{
  { 
    # stopCluster(cl)
    cl <- makeCluster(detectCores() - 1)
    clusterEvalQ(cl, {library(tidyverse)})
  }
  
  temp <- read.vcfR(file = "/home/jinoo/ExAC.r0.3.1.sites.vep.vcf", verbose = T) 
  temp_tidy <- temp@fix %>% as_tibble() %>% select(-ID, -QUAL, -FILTER, -REF, -ALT)
  AC_AF <- parLapply(cl = cl, X = temp_tidy$INFO, fun = function(x){
    temp_info <- x %>% str_split(string = ., pattern = ";", simplify = T) %>% 
      str_split(string = ., pattern = "=", simplify = T) %>%
      as_tibble() %>% mutate(V2 = as.numeric(V2)) %>% 
      subset.data.frame(., subset = (V1 == "AC_Adj" | V1 == "AN_Adj"))
    
    return(tibble(AF_adj = temp_info$V2[1]/temp_info$V2[2]))
  }) %>% bind_rows()
  
  AC_AF <- bind_rows(AC_AF)
  temp_tidy <- select(temp_tidy, -INFO)
  ExAC_MAF <- bind_cols(temp_tidy, AC_AF)
  
  fwrite(x = ExAC_MAF, file = "../ExAC_MAF.txt",sep = "\t", row.names = F, col.names = T, quote = F)
}

  

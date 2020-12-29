excel_go_geneset <- function(geneset){
  max_len <- geneset[[1]] %>% lapply(X = ., function(x){
    length(x) %>% return()
  }) %>% unlist() %>% max()
  
  DF <- geneset[[1]] %>% lapply(X = ., function(x){
    temp <- vector(mode = "character", max_len)
    for(index in 1:(length(x))){
      if(!is.na(x[index]))
        temp[index] <- x[index]
    }
    
    temp %>% as_tibble() %>% return()
  }) %>% bind_cols()
  colnames(DF) <- geneset[[2]]
  
  return(DF)
}
excel_go_geneset(geneset_all)


temp_list <- list()
temp_list[[1]] <- fread("WES_merge_SKATO_GO_set2.txt") %>% 
  separate(SetID_003, into = c("types_of_variants", "KEY"), sep = "__") %>% 
  select(-SetID_001)
temp_list[[2]] <- fread("NeuroX_SKATO_GO_set2.txt") %>%
  separate(SetID_003, into = c("types_of_variants", "KEY"), sep = "__") %>% 
  select(-SetID_001)
temp_list[[3]] <- fread("PPMI_SKATO_GO_set2.txt") %>%
  separate(SetID_003, into = c("types_of_variants", "KEY"), sep = "__") %>% 
  select(-SetID_001)

for(value in geneset[[2]]){
  a <- temp_list[[1]] %>% filter(KEY == value) %>% bind_cols(tibble(Data = "NINDS(exome)"), .) %>% 
    select(-KEY)
  b <- temp_list[[2]] %>% filter(KEY == value) %>% bind_cols(tibble(Data = "NeuroX"), .) %>% 
    select(-KEY)
  c <- temp_list[[3]] %>% filter(KEY == value) %>% bind_cols(tibble(Data = "PPMI"), .) %>% 
    select(-KEY)
  
  write_delim(x = value %>% as_tibble(), path = "temp.txt", delim = "\t", append = T)
  write_delim(x = bind_rows(a,b,c), path = "temp.txt", delim = "\t", append = T)
  write_delim(x = tibble(blank = " "), path = "temp.txt", delim = "\t", append = T)
}





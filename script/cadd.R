temp <- fread(file = "../all_case_106.vcf", skip = 133, header = T) %>%
  as_tibble() %>% select(1:8)
dir.create(path = "./cadd_vcf")
setwd("cadd_vcf/")
n <- 1
last_len <- nrow(temp)
while(T){
  if(n > last_len)
    break
  temp %>% slice(n:(n+40000)) %>%
    write_delim(path = paste0("cadd_vcf_", n, ".vcf"), delim = "\t", col_names = T)
  n <- n + 40000
}

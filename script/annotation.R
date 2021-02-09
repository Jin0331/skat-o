## annotation

library(tidyverse)
library(data.table)
library(parallel)
library(doMC)
library(vcfR)
library(glue)

setwd("data/anno/")

### vcf split
### exclude_case 138
###  all_case 133
core <- 7
temp_tidy <- fread(input = "../all_case_106.vcf", sep = "\t", quote = "\n", skip = 133, header = T) %>% as_tibble()

vcf_list <- list()
vcf_list[1] <- list(temp_tidy[1:round(nrow(temp_tidy)/core),])
for(i in 2:(core)){
  if(i == core){
    vcf_list[i] <- list(temp_tidy[((nrow(vcf_list[[1]])*(i-1))+1):(nrow(temp_tidy)),])
    gc()
    break;
  }
  vcf_list[i] <- list(temp_tidy[((nrow(vcf_list[[1]])*(i-1))+1):(nrow(vcf_list[[1]])*i),])
  gc()
}
rm(temp_tidy);gc()


## for annotation, vcf out
for(i in 1:length(vcf_list)){
  fwrite(x = vcf_list[[i]], file = paste0("vcf",i,".vcf"), row.names = F, col.names = T, sep = "\t", quote = F, nThread = 14)
  print(i)
};rm(vcf_list);gc()

## vcf annotation
path <- list.files(full.names = T)
registerDoMC(core)
system.time(
  foreach(i=1:core) %dopar% {
    system(glue("/tools/annovar/table_annovar.pl {path[i]} /tools/annovar/humandb/humandb -buildver hg19 -out {path[i]} -remove -protocol knownGene,avsnp150 -operation g,f -nastring . -vcfinput"))
  }
);gc() #

#### make annotated fix
system("rm -rf *.txt");system("rm -rf *.avinput")
temp <- list.files(pattern = "hg19")

mul_vcf <- mclapply(temp, read.vcfR, convertNA =T, checkFile = F, mc.cores = 2)
mul_anno_fix <- lapply(X = mul_vcf, FUN = function(temp){
  temp2 <- vcfR2tidy(temp, info_only = T, single_frame = F, toss_INFO_column = T)
  return(temp2$fix)
})
mul_anno_fix <- bind_rows(mul_anno_fix)
rm(mul_vcf)
fwrite(x = mul_anno_fix, file = "../all_case_106_fix_raw.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n")

###  CADD annotation
file_path <- list.files()

temp <- lapply(X = file_path, function(value){
  fread(value, header = T) %>% as_tibble() %>% 
    mutate(`#CHROM` = ifelse(`#CHROM` == "X", "23", ifelse(`#CHROM` == "Y", "24", `#CHROM`))) %>%
    mutate(`#CHROM` = as.integer(`#CHROM`)) %>% 
    rename(CHROM = `#CHROM`)
}) %>% bind_rows() %>% unique()

gnomad <- fread("https://s3.us-west-2.amazonaws.com/secure.notion-static.com/bb4a5528-57ef-4500-bca4-c26518c0329a/GnomAD_V2.txt?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAT73L2G45O3KS52Y5%2F20210208%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20210208T064830Z&X-Amz-Expires=86400&X-Amz-Signature=3181a84c5f781a52d664c42c31bcd1dbfd29e59cd576ffee5dce6dbef123e751&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22GnomAD_V2.txt%22",
                header = T) %>% as_tibble() %>% 
  select(-FILTER)

raw_fix <- fread("../exclude_case_fix_raw.txt") %>% as_tibble() %>% 
  mutate(`CHROM` = ifelse(`CHROM` == "X", "23", 
                           ifelse(`CHROM` == "Y", "24", `CHROM`))) %>%
  mutate(`CHROM` = as.integer(`CHROM`)) %>% 
  left_join(x = ., y = gnomad, by = c("CHROM", "POS", "REF", "ALT")) %>% 
  left_join(x = ., y = temp, by = c("CHROM", "POS", "REF", "ALT")) %>% 
  mutate(AF = ifelse(is.na(AF), 0.009, AF), AF_eas = ifelse(is.na(AF_eas), 0.009, AF_eas), 
         AF_eas_kor = ifelse(is.na(AF_eas_kor), 0.009, AF_eas_kor))

raw_fix %>% write_delim(file = "../exclude_case_fix_anno.txt", col_names = T, delim = "\t")



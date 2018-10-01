install.packages("glue")
install.packages("vcfR")
install.packages("data.table")
install.packages("foreach")
install.packages("doMC")
install.packages("dplyr")
library(glue);library(vcfR);library(data.table);library(foreach);library(doMC)
library(dplyr)

################################# row data #################################
### vcf-merge
## before vcftool using
system("sed -i 's/chr//g' /home/jinoo/skat-o/row_data/NABEC_EXOME_cohort_control.vcf");system("sed -i 's/chr//g' /home/jinoo/skat-o/row_data/coriell_pd_exomes_case.vcf")
system("bgzip /home/jinoo/skat-o/row_data/NABEC_EXOME_cohort_control.vcf");system("bgzip /home/jinoo/skat-o/row_data/coriell_pd_exomes_case.vcf")
system("tabix -p vcf /home/jinoo/skat-o/row_data/NABEC_EXOME_cohort_control.vcf.gz");system("tabix -p vcf /home/jinoo/skat-o/row_data/coriell_pd_exomes_case.vcf.gz")
system("vcf-merge --remove-duplicates /home/jinoo/skat-o/row_data/NABEC_EXOME_cohort_control.vcf.gz /home/jinoo/skat-o/row_data/coriell_pd_exomes_case.vcf.gz | bgzip -c > /home/jinoo/skat-o/row_data/merge_row.vcf.gz")

test <- read.vcfR(file = "merge_row.vcf.gz", convertNA = T, checkFile = F)
test <- vcfR2tidy(test, info_only = T, single_frame = F, toss_INFO_column = T)
test_fix <- test$fix;rm(test)
fwrite(x = test_fix, file = "/home/skat-0/row_data/vcf_bind_0930_row_fix.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n");rm(test_fix) ## fix out
################################# row data #################################



#### directory make
setwd("~/skat-o/")
date <- "1001"; date <- paste0(date, "_test") ######## date set, input !!!!!!!
system(glue("mkdir {test_date}",test_date = date)) 
setwd(paste0(getwd(), "/",date))

#### vcf QC -- KGGSeq
system("java -Xmx10g -jar /home/lee/kggseq10hg19/kggseq.jar --vcf-file /home/jinoo/skat-o/row_data/merge_row.vcf.gz --ped-file /home/jinoo/skat-o/skato_0918.ped --out skatQC_KGGSeq --o-vcf --gty-qual 20 --gty-dp 8 --gty-sec-pl 20 --vcf-filter-in PASS --seq-qual 50 --seq-mg 20 --seq-fs 60 --hwe-control 1E-5 --nt 7")
system("gzip -d skatQC_KGGSeq.flt.vcf.gz") ## bgzip and tabix
system("bgzip skatQC_KGGSeq.flt.vcf")
system("tabix -p vcf skatQC_KGGSeq.flt.vcf.gz")
test <- read.vcfR(file = "skatQC_KGGSeq.flt.vcf.gz", convertNA = T, checkFile = F)
test <- vcfR2tidy(test, info_only = T, single_frame = F, toss_INFO_column = T)
test_fix <- test$fix;rm(test)
fwrite(x = test_fix, file = "skatQC_KGGSeq_GVQC.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n");rm(test_fix) ## fix out
rm(test_fix)

#### vcf QC -- vcftools
system("vcftools --gzvcf skatQC_KGGSeq.flt.vcf.gz --min-alleles 2 --max-alleles 2 --max-missing 0.8 --recode --out QC_80_allele_individual_vcftools_flt")
test <- read.vcfR(file = "QC_80_allele_individual_vcftools_flt", convertNA = T, checkFile = F)
test <- vcfR2tidy(test, info_only = T, single_frame = F, toss_INFO_column = T)
test_fix <- test$fix;rm(test)
fwrite(x = test_fix, file = "skatQC_allele_individual_flt_vcftools.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n");rm(test_fix) ## fix out


## row merge & recode vcf load(7 div)
temp <- fread(input = "QC_80_allele_individual_vcftools_flt.recode.vcf", sep = "\t", quote = "\n", skip = 210, header = T, nThread = 14)
temp_tidy <- dplyr::tbl_df(temp);rm("temp")
vcf_list <- list()
vcf_list[1] <- list(temp_tidy[1:round(nrow(temp_tidy)/core),])
for(i in 2:(core)){
  if(i == core){
    vcf_list[i] <- list(temp_tidy[((nrow(vcf1)*(i-1))+1):(nrow(temp_tidy)),])
    break;
  }
  vcf_list[i] <- list(temp_tidy[((nrow(vcf1)*(i-1))+1):(nrow(vcf1)*i),])
}
rm(temp_tidy);rm(vcf1)

## for annotation, vcf out
system("mkdir vcf");setwd("vcf/")
for(i in 1:length(vcf_list)){
  fwrite(x = vcf_list[[i]], file = paste0("vcf",i,".vcf"), row.names = F, col.names = T, sep = "\t", quote = F, nThread = 14)
  print(i)
}
rm(vcf_list);gc()

## vcf annotation
path <- list.files(full.names = T)
registerDoMC(core)
system.time(
  foreach(i=1:7) %dopar% {
    system(glue("/home/lee/annovar/table_annovar.pl {path[i]} /home/lee/annovar/humandb/ -buildver hg19 -out {path[i]} -remove -protocol knownGene,refGene,avsnp147,cadd13gt10 -operation g,g,f,f -nastring . -vcfinput"))
  }
);gc() ##7000_ 3div

temp <- list.files(pattern = "*hg19")
mul_vcf <- mclapply(temp, read.vcfR, convertNA =T, checkFile = F, mc.cores = 2)

vcf_bind1 <- rbind2(vcf_1,vcf_2);rm(vcf_1,vcf_2)
vcf_bind2 <- rbind2(vcf_3,vcf_4);rm(vcf_3,vcf_4)
vcf_bind3 <- rbind2(vcf_bind1,vcf_bind2);rm(vcf_bind1,vcf_bind2)
vcf_bind4 <- rbind2(vcf_5,vcf_6);rm(vcf_5,vcf_6)
vcf_bind5 <- rbind2(vcf_bind3,vcf_bind4);rm(vcf_bind3,vcf_bind4)
vcf_bind_complete <- rbind2(vcf_bind5, vcf_7);rm(vcf_7,vcf_bind5)

# annotation fix & vcf out
test <- vcfR2tidy(vcf_bind_complete, info_only = T, single_frame = F, toss_INFO_column = T)
test_fix <- test$fix;rm(test)
fwrite(x = test_fix, file = "annotation_fix.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n") ## fix out
write.vcf(x = vcf_bind_complete, file = "all_c_QC_80.vcf.gz",APPEND = F);rm(vcf_bind_complete) ## vcf out
system("gzip -d all_c_QC_80.vcf.gz") ## bgzip and tabix
system("bgzip all_c_QC_80.vcf")
system("tabix -p vcf all_c_QC_80.vcf.gz")
system("rm -rf ")



################### skat-o preprocessing #############################
test_fix <- read.table(file = "annotation_fix.txt", sep = "\t", header = T, stringsAsFactors = F)
geneset <- read.csv("/home/jinoo/skat-o/parkinson_genset.txt", stringsAsFactors = F,header = F)
geneset <- as.character(geneset[,1])

## NULL dataframe
variant_all_parkinson <- data.frame(CHROM = NA, POS = NA, ID = NA, dbSNP147 = NA, REF = NA , ALT = NA, Gene.knownGene = NA, ExonicFunc.knownGene = NA, CADD13_PHRED=NA)[numeric(0), ]

## Parkinson geneset variant
# 1:17312586_G/A   grp format

for(i in 1:length(geneset)){
  variant_all_parkinson <- rbind(variant_all_parkinson, subset(test_fix, subset = ((test_fix$Gene.knownGene == geneset[i] & test_fix$Func.knownGene == "exonic")), 
                                                               select = c("CHROM", "POS", "ID", "avsnp147","REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED")), fill=T)
}

### 1. nonsynonymous geneset
nonsynonymous <- subset(variant_all_parkinson, subset = ((variant_all_parkinson$ExonicFunc.knownGene == "nonsynonymous_SNV")) | variant_all_parkinson$ExonicFunc.knownGene == "stopgain" | variant_all_parkinson$ExonicFunc.knownGene ==  "stoploss" 
                            | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_deletion" | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_insertion" 
                            | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_block_substitution" | variant_all_parkinson$ExonicFunc.knownGene ==  "splicing")

i <- NULL;nonsynonymous_geneset <- c()
for(i in 1:nrow(nonsynonymous)){
  nonsynonymous_geneset <- c(nonsynonymous_geneset, paste0(nonsynonymous[i,1], ":", nonsynonymous[i,2], "_", nonsynonymous[i,5], "/", nonsynonymous[i,6]))
}

### 2. CADD > 12.37 variant
cadd <- subset(variant_all_parkinson, subset = ( ((variant_all_parkinson$ExonicFunc.knownGene == "nonsynonymous_SNV")) & variant_all_parkinson$CADD13_PHRED > 12.37))
i <- NULL;cadd_geneset <- c()
for(i in 1:nrow(cadd)){
  cadd_geneset <- c(cadd_geneset, paste0(cadd[i,1], ":", cadd[i,2], "_", cadd[i,5], "/", cadd[i,6]))
}

### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
lof <- subset(variant_all_parkinson, subset = ( variant_all_parkinson$ExonicFunc.knownGene == "stopgain" | variant_all_parkinson$ExonicFunc.knownGene ==  "stoploss" 
                                      | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_deletion" | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_insertion" 
                                      | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_block_substitution" | variant_all_parkinson$ExonicFunc.knownGene ==  "splicing"))
i <- NULL;lof_geneset <- c()
for(i in 1:nrow(lof)){
  lof_geneset <- c(lof_geneset, paste0(lof[i,1], ":", lof[i,2], "_", lof[i,5], "/", lof[i,6]))
}

nonsynonymous_geneset <- t(data.frame(nonsynonymous_geneset, stringsAsFactors = F));rownames(nonsynonymous_geneset) <- "nonsynonymous"
cadd_geneset <- t(data.frame(cadd_geneset, stringsAsFactors = F));rownames(cadd_geneset) <- "cadd"
lof_geneset <- t(data.frame(lof_geneset, stringsAsFactors = F));rownames(lof_geneset) <- "lof"
write.table(x = nonsynonymous_geneset, "skat_grp.grp", sep = "\t",row.names = T, quote = F, col.names = F, append = T)
write.table(x = cadd_geneset, "skat_grp.grp", sep = "\t", row.names = T, quote = F, col.names = F, append = T)
write.table(x = lof_geneset, "skat_grp.grp", sep = "\t", row.names = T, quote = F, col.names = F, append = T)


#### gene subset
### exonic
variant_all_parkinson_gene <- unique(c(variant_all_parkinson$Gene.knownGene, variant_all_parkinson$Gene.refGene))
paste_temp <- NULL
temp <- NULL
for( i in 1:length(variant_all_parkinson_gene)){
  paste_temp <- NULL
  temp <- subset(variant_all_parkinson, subset = ((Gene.knownGene %in% variant_all_parkinson_gene[i])))
  temp <- subset(temp, subset = ((temp$ExonicFunc.knownGene == "nonsynonymous_SNV"))| temp$ExonicFunc.knownGene == "stopgain"
             | temp$ExonicFunc.knownGene ==  "stoploss" | temp$ExonicFunc.knownGene ==  "frameshift_deletion"
             | temp$ExonicFunc.knownGene ==  "frameshift_insertion" | temp$ExonicFunc.knownGene ==  "frameshift_block_substitution"
             | temp$ExonicFunc.knownGene ==  "splicing")
  if(nrow(temp) != 0){
  for(j in 1:nrow(temp)){
    paste_temp <- c(paste_temp, paste0(temp[j,1], ":", temp[j,2],"_", temp[j,5],"/",temp[j,6]))
  }
  paste_temp <- t(data.frame(paste_temp,stringsAsFactors = F));rownames(paste_temp) <- variant_all_parkinson_gene[i]
  write.table(paste_temp, "variant_all_parkinson_gene_geneset_gene.grp", sep = "\t", row.names = T, col.names = F,append = T, quote = F)
  }
}



### skat-ot test
system("pwd")
system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf skat_grp.grp  -out test_0928_maf005.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.05 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2")
system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf skat_grp.grp  -out test_0928_maf003.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.03 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2")
system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf skat_grp.grp  -out test_0928_maf001.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.01 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2")

system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf variant_all_parkinson_gene_geneset_gene.grp  -out test_0928_exonic_all_gene005.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.05 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2") 
system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf variant_all_parkinson_gene_geneset_gene.grp  -out test_0928_exonic_all_gene003.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.03 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2") 
system("/home/lee/epacts_0913/bin/epacts-group -vcf vcf_bind_0917_gz.vcf.gz -groupf variant_all_parkinson_gene_geneset_gene.grp  -out test_0928_exonic_all_gene001.skat -ped /home/jinoo/skat-o/skato_0918.ped -max-maf 0.01 -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2")


### result_adjusted_p value
file <- list.files(pattern = "test*")
for(i in 1:length(file)){
  temp <- read.table(file = file[i], header = F)
  system(glue("rm -rf {remove}", remove = file[i]))
  temp_colname <- c("#CHROM","BEGIN","END","MARKER_ID","NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE","STATRHO");colnames(temp) <- temp_colname
  PVALUE_ADJUSTED <- p.adjust(temp$PVALUE, p.adjust.methods[4])
  temp <- cbind(temp, PVALUE_ADJUSTED)
  write.table(x = temp, file = file[i], col.names = T, row.names = F, quote = F, sep = "\t")
}

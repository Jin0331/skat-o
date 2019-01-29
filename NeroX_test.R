library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(dplyr)
# tool path 
Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/lee/tool/:/home/lee/annovar/:/home/lee/epacts_0913/bin/")
system("")
## unknown snp QC(e.g. minor | major 0)
bim_qc <- read.table(file = "dbGaP_NeuroX_filter.bim", header = F, stringsAsFactors = F)
bim_zero_snp <- subset.data.frame(x = bim_qc, subset = (V5 == '0' | V6 =='0'), select = "V2")[,1]
write.table(x = bim_zero_snp, file = "remove_snp_file_0129.txt", sep = "\t", quote = F, row.names = F, col.names = F)

system("mkdir data");setwd("data/")
system("plink --bfile ../dbGaP_NeuroX_filter --exclude ../remove_snp_file_0129.txt --make-bed --out skatQC")
system("./update_build.sh skatQC ../NeuroX_15036164_A-b37.Ilmn.strand skatQC_strand") ### strand update
system("plink --bfile skatQC_strand --reference-allele ../NeuroX_15036164_A-b37.strand.RefAlt --make-bed --out skatQC_strand_ref")
setwd("../")

## sample QC
system("mkdir plink");setwd("plink/")
system("plink --bfile ../data/skatQC_strand_ref --het --out skatQC_het") ### heterogygosity
system("plink --bfile ../data/skatQC_strand_ref --indep-pairwise 50 5 0.2 --out skat_pruning")
system("plink --bfile ../data/skatQC_strand_ref --genome full --out skatQC_IBD") ### MDS
system("plink --bfile ../data/skatQC_strand_ref --read-genome skatQC_IBD.genome --extract skat_pruning.prune.in --mds-plot 20 --cluster --out skatQC_MDS")
system("plink --bfile ../data/skatQC_strand_ref --missing --out skatQC_missing") ### individual missing
system("plink --bfile ../data/skatQC_strand_ref --check-sex --out skatQC_sex")

#1 heterozygosity
plink_het <- read.table(file = "skatQC_het.het",header = T, stringsAsFactors = F)
F_mean_4SD <- mean(plink_het$F) - (sd(plink_het$F)*4)
removal_het <- subset(plink_het, subset = (plink_het$F < F_mean_4SD), select = "IID")[,1] ### het removal check)

#2 population outliers - MDS
plink_MDS <- read.table(file = "skatQC_MDS.mds", header = T, stringsAsFactors = F)
col_mean <- apply(plink_MDS[,4:23], 2, mean)
col_sd <- apply(plink_MDS[,4:23], 2, sd)
removal_MDS <- NULL
for(i in 1:nrow(plink_MDS)){
  for( j in 4:23){ 
    if( plink_MDS[i,j] > (col_mean[j-3] + 4*col_sd[j-3]) | plink_MDS[i,j] < (col_mean[j-3] - 4*col_sd[j-3]) )
      removal_MDS <- c(removal_MDS, plink_MDS[i,1])
  }
}

#3 missing
plink_missing <- read.table(file = "skatQC_missing.imiss", header = T, stringsAsFactors = F)
removal_imissing <- subset(plink_missing, subset = (F_MISS > 0.15), select = c("IID","F_MISS"))[,1]

#4 gender
plink_gender <- read.table(file = "skatQC_sex.sexcheck", header = T, stringsAsFactors = F)
removal_gender <- subset(plink_gender, subset = (STATUS == "PROBLEM"), select = c("IID"))[,1]

#5 vcf QC -- vcftools, ajk value

# system("pseq NeuroX new-project --resources /home/lee/pseq/hg19/")
# system("pseq NeuroX.pseq load-plink --file skatQC_snp --id skatQC_snp")
# system("pseq NeuroX.pseq write-vcf > NeuroX_filter.vcf")
# system("pseq NeuroX_filter.vcf write-vcf --format BGZF --file NeuroX_filter.vcf.gz")
# system("mv NeuroX_filter.vcf NeuroX_filter.vcf.bak")
# system("gzip -d NeuroX_filter.vcf.gz")

system("plink --bfile ../data/skatQC_strand_ref --recode vcf-iid --allow-no-sex --keep-allele-order --out NeuroX_filter")
system("sed -i 's/chr//g' NeuroX_filter.vcf")

#6 merge & write
removal_iid <- unique(c(removal_MDS,removal_het,removal_imissing, removal_gender)) #, removal_ajk
write.table(x = removal_iid, file = "plink_ind.txt",sep = "\n", quote = F, row.names = F, col.names = F)
rm(list=ls());gc()

system("plink --bfile ../data/skatQC_strand_ref --geno 0.01 --hwe 1e-6 --make-bed --out skatQC_strand_ref_snp --noweb")
system("plink --bfile skatQC_strand_ref_snp --recode vcf-iid --allow-no-sex --keep-allele-order --out ../NeuroX_filter_qc")
system("sed -i 's/chr//g' NeuroX_filter_qc.vcf")
system("bgzip NeuroX_filter_qc.vcf")
system("tabix -p vcf NeuroX_filter_qc.vcf.gz") 

system("vcftools --gzvcf NeuroX_filter_qc.vcf.gz --min-alleles 2 --max-alleles 2 --remove plink/plink_ind.txt --recode --out NeuroX_vcftool_qc.flt")

### vcf split
core <- 9
temp <- fread(input = "NeuroX_vcftool_qc.flt.recode.vcf", sep = "\t", quote = "\n", skip = 30, header = T, nThread = 5) ### if qc, skip = 3 and select remove

temp_tidy <- dplyr::tbl_df(temp)
rm("temp")
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
system("mkdir vcf")
setwd("vcf/")
for(i in 1:length(vcf_list)){
  fwrite(x = vcf_list[[i]], file = paste0("vcf",i,".vcf"), row.names = F, col.names = T, sep = "\t", quote = F, nThread = 14)
  print(i)
}
rm(vcf_list);gc()


## vcf annotation
path <- list.files(full.names = T)
registerDoMC(9)
system.time(
  foreach(i=1:core) %dopar% {
    system(glue("table_annovar.pl {path[i]} /home/lee/annovar/humandb/ -buildver hg19 -out {path[i]} -remove -protocol knownGene,avsnp147,cadd13gt10 -operation g,f,f -nastring . -vcfinput"))
  }
);gc() #

#### make annotated fix
system("rm -rf *.txt");system("rm -rf *.avinput")
temp <- list.files(pattern = "*hg19")

mul_vcf <- mclapply(temp, read.vcfR, convertNA =T, checkFile = F, mc.cores = 2)
mul_anno_fix <- lapply(X = mul_vcf, FUN = function(temp){
  temp2 <- vcfR2tidy(temp, info_only = T, single_frame = F, toss_INFO_column = T)
  return(temp2$fix)
})
mul_anno_fix <- bind_rows(mul_anno_fix)
rm(mul_vcf)
fwrite(x = mul_anno_fix, file = "../NeuroX_annotation_fix_0129.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n")

setwd("../")
system("bgzip NeuroX_vcftool_qc.flt.recode.vcf")
system("tabix -p vcf NeuroX_vcftool_qc.flt.recode.vcf.gz")
################### skat-o preprocessing #############################
######################################################################



## new ped
setwd("../")
skat_row_ped <- read.table(file = "../NeuroX_filter_epact_row.ped", header = T, stringsAsFactors = F)
skat_subset <- read.table(file = "plink_ind.txt", header = F, stringsAsFactors = F)[,1]
'%!in%' <- function(x,y)!('%in%'(x,y))
skat_QC_ped <- NULL
skat_QC_ped <- subset(skat_row_ped, subset = skat_row_ped$fid %!in% skat_subset)
write.table(x = skat_QC_ped, file = "../NeuroX_epacts_0129.ped", col.names = T, row.names = F, sep = "\t", quote = F)

system("mkdir result");setwd("result")
system("pwd")
test_fix <- read.table(file = "../NeuroX_annotation_fix_0129.txt", sep = "\t", header = T, stringsAsFactors = F)
geneset <- read.csv("/home/jinoo/skat-o/parkinson_genset.txt", stringsAsFactors = F,header = F)
# geneset <- read.csv("/home/jinoo/skat-o/LSD_geneset.txt", stringsAsFactors = F,header = F)
geneset <- as.character(geneset[,1])

## NULL dataframe
variant_all_parkinson <- data.frame(CHROM = NA, POS = NA, ID = NA, snp131 = NA, REF = NA , ALT = NA, Gene.knownGene = NA, ExonicFunc.knownGene = NA, CADD13_PHRED=NA)[numeric(0), ]

## Parkinson geneset variant
# 1:17312586_G/A   grp format

for(i in 1:length(geneset)){
  variant_all_parkinson <- rbind(variant_all_parkinson, subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")), 
                                                               select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED")))
}

#### lapply convert !!!!!!!

### for gene, variant number
variant_num1 <- NULL;variant_num2 <- NULL;data_temp <- NULL; gene_to_variant <- NULL;
gene_to_variant <- data.frame(Gene = NA, variant_count=NA)[numeric(0), ]
for(i in 1:length(geneset)){
  variant_num1 <- subset(variant_all_parkinson, subset = ( (Gene.knownGene == geneset[i] & CADD13_PHRED > 12.37))) ### for geneset , all-variant
  variant_num2 <- nrow(subset(variant_num1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" | ExonicFunc.knownGene == "stopgain"
                                                      | ExonicFunc.knownGene ==  "stoploss" | ExonicFunc.knownGene ==  "frameshift_deletion" 
                                                      | ExonicFunc.knownGene ==  "frameshift_insertion" 
                                                      | ExonicFunc.knownGene ==  "frameshift_block_substitution" | ExonicFunc.knownGene ==  "splicing")))
  data_temp <- data.frame(Gene = geneset[i], variant_num2)
  gene_to_variant <- rbind(gene_to_variant, data_temp) #### result
}

View(gene_to_variant)
write.table(gene_to_variant, file = "../gene-variant_count_parkinson_0128.txt", quote = F, sep = "\t")

### 1. nonsynonymous geneset
nonsynonymous <- subset(variant_all_parkinson, subset = ((variant_all_parkinson$ExonicFunc.knownGene == "nonsynonymous_SNV")) | variant_all_parkinson$ExonicFunc.knownGene == "stopgain" | variant_all_parkinson$ExonicFunc.knownGene ==  "stoploss" 
                        | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_deletion" | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_insertion" 
                        | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_block_substitution" | variant_all_parkinson$ExonicFunc.knownGene ==  "splicing")

i <- NULL;nonsynonymous_geneset <- c()
for(i in 1:nrow(nonsynonymous)){
  nonsynonymous_geneset <- c(nonsynonymous_geneset, paste0(nonsynonymous[i,1], ":", nonsynonymous[i,2], "_", nonsynonymous[i,4], "/", nonsynonymous[i,5]))
}

### 2. CADD > 12.37 variant
cadd <- subset(variant_all_parkinson, subset = ( ((variant_all_parkinson$ExonicFunc.knownGene == "nonsynonymous_SNV")) & variant_all_parkinson$CADD13_PHRED > 12.37))
i <- NULL;cadd_geneset <- c()
for(i in 1:nrow(cadd)){
  cadd_geneset <- c(cadd_geneset, paste0(cadd[i,1], ":", cadd[i,2], "_", cadd[i,4], "/", cadd[i,5]))
}

### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
lof <- subset(variant_all_parkinson, subset = ( variant_all_parkinson$ExonicFunc.knownGene == "stopgain" | variant_all_parkinson$ExonicFunc.knownGene ==  "stoploss" 
                                                | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_deletion" | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_insertion" 
                                                | variant_all_parkinson$ExonicFunc.knownGene ==  "frameshift_block_substitution" | variant_all_parkinson$ExonicFunc.knownGene ==  "splicing") & variant_all_parkinson$CADD13_PHRED > 12.37)
i <- NULL;lof_geneset <- c()
if(nrow(lof) != 0 ){
  for(i in 1:nrow(lof)){
    lof_geneset <- c(lof_geneset, paste0(lof[i,1], ":", lof[i,2], "_", lof[i,4], "/", lof[i,5]))
  }
}
nonsynonymous_geneset <- t(data.frame(nonsynonymous_geneset, stringsAsFactors = F));rownames(nonsynonymous_geneset) <- "nonsynonymous"
cadd_geneset <- t(data.frame(cadd_geneset, stringsAsFactors = F));rownames(cadd_geneset) <- "cadd"
if(nrow(lof) !=0){
    lof_geneset <- t(data.frame(lof_geneset, stringsAsFactors = F));rownames(lof_geneset) <- "lof"}
system("rm -rf skat_grp.grp")
write.table(x = nonsynonymous_geneset, "skat_grp.grp", sep = "\t",row.names = T, quote = F, col.names = F, append = T)
write.table(x = cadd_geneset, "skat_grp.grp", sep = "\t", row.names = T, quote = F, col.names = F, append = T)
write.table(x = lof_geneset, "skat_grp.grp", sep = "\t", row.names = T, quote = F, col.names = F, append = T)


#### gene subset
### exonic
variant_all_parkinson_gene <- unique(c(variant_all_parkinson$Gene.knownGene))
paste_temp <- NULL
temp <- NULL
system("rm -rf variant_all_parkinson_gene_geneset_gene.grp")
for( i in 1:length(variant_all_parkinson_gene)){
  paste_temp <- NULL
  temp <- subset(variant_all_parkinson, subset = ((Gene.knownGene %in% variant_all_parkinson_gene[i])) & CADD13_PHRED > 12.37)
  temp <- subset(temp, subset = ((temp$ExonicFunc.knownGene == "nonsynonymous_SNV"))| temp$ExonicFunc.knownGene == "stopgain"
                 | temp$ExonicFunc.knownGene ==  "stoploss" | temp$ExonicFunc.knownGene ==  "frameshift_deletion"
                 | temp$ExonicFunc.knownGene ==  "frameshift_insertion" | temp$ExonicFunc.knownGene ==  "frameshift_block_substitution"
                 | temp$ExonicFunc.knownGene ==  "splicing")
  if(nrow(temp) != 0){
    for(j in 1:nrow(temp)){
      paste_temp <- c(paste_temp, paste0(temp[j,1], ":", temp[j,2],"_", temp[j,4],"/",temp[j,5]))
    }
    paste_temp <- t(data.frame(paste_temp,stringsAsFactors = F));rownames(paste_temp) <- variant_all_parkinson_gene[i]
    write.table(paste_temp, "variant_all_parkinson_gene_geneset_gene.grp", sep = "\t", row.names = T, col.names = F,append = T, quote = F)
  }
}



### skat-ot test

MAF <- c(0.01, 0.03)
for(i in MAF){
  system(glue("epacts-group -vcf ../NeuroX_vcftool_qc.flt.recode.vcf.gz -groupf skat_grp.grp -out NeuroX_0129_maf{maf}_PARK.skat -ped ../NeuroX_epacts_0129.ped -max-maf {maf} -pheno disease -cov sex -missing ./. -test skat -skat-o -run 4", maf = i))
  system(glue("epacts-group -vcf ../NeuroX_vcftool_qc.flt.recode.vcf.gz -groupf variant_all_parkinson_gene_geneset_gene.grp -out NeuroX_0129_maf{maf}_PARK_gene.skat -ped ../NeuroX_epacts_0129.ped -max-maf {maf} -pheno disease -cov sex -missing ./. -test skat -skat-o -run 4", maf = i))
}

### result_adjusted_p value
system("find . ! -name '*.epacts' -delete")
file <- list.files()

for(i in 1:length(file)){
  temp <- read.table(file = file[i], header = F)
  system(glue("rm -rf {remove}", remove = file[i]))
  temp_colname <- c("#CHROM","BEGIN","END","MARKER_ID","NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE","STATRHO");colnames(temp) <- temp_colname
  PVALUE_ADJUSTED <- p.adjust(temp$PVALUE, p.adjust.methods[4])
  temp <- cbind(temp, PVALUE_ADJUSTED)
  write.table(x = temp, file = file[i], col.names = T, row.names = F, quote = F, sep = "\t")
}

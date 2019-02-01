install.packages("SKAT")
library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(dplyr)
library(SKAT)

# tool path 
Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/lee/tool/:/home/lee/annovar/:/home/lee/epacts_0913/bin/")

## unknown snp QC(e.g. minor | major 0)
bim_qc <- read.table(file = "dbGaP_NeuroX_filter.bim", header = F, stringsAsFactors = F)
bim_zero_snp <- subset.data.frame(x = bim_qc, subset = (V5 == '0' | V6 =='0'), select = "V2")[,1]
write.table(x = bim_zero_snp, file = "remove_snp_file_0129.txt", sep = "\t", quote = F, row.names = F, col.names = F)

system("mkdir plink");setwd("plink/")
system("plink --bfile ../dbGaP_NeuroX_filter --exclude ../remove_snp_file_0129.txt --make-bed --out skatQC")
system("../update_build.sh skatQC ../NeuroX_15036164_A-b37.Ilmn.strand skatQC_strand") ### strand update
system("plink --bfile skatQC_strand --reference-allele ../NeuroX_15036164_A-b37.strand.RefAlt --make-bed --out skatQC_strand_ref")
setwd("../")

## sample QC
system("plink --bfile skatQC_strand_ref --het --out skatQC_het") ### heterogygosity
system("plink --bfile skatQC_strand_ref --indep-pairwise 50 5 0.2 --out skat_pruning")
system("plink --bfile skatQC_strand_ref --genome full --out skatQC_IBD") ### MDS
system("plink --bfile skatQC_strand_ref --read-genome skatQC_IBD.genome --extract skat_pruning.prune.in --mds-plot 20 --cluster --out skatQC_MDS")
system("plink --bfile skatQC_strand_ref --missing --out skatQC_missing") ### individual missing
system("plink --bfile skatQC_strand_ref --check-sex --out skatQC_sex")

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
system("plink --bfile skatQC_strand_ref --recode vcf-iid --allow-no-sex --keep-allele-order --out NeuroX_filter")
system("sed -i 's/chr//g' NeuroX_filter.vcf")

#6 merge & write
removal_iid <- unique(c(removal_MDS,removal_het,removal_imissing, removal_gender)) #, removal_ajk
removal_iid <- data.frame(FID=removal_iid, IID=removal_iid, stringsAsFactors = F)
write.table(x = removal_iid, file = "plink_ind.txt",sep = "\t", quote = F, row.names = F, col.names = F)
rm(list=ls());gc()

#7 snp qc and write vcf format
system("plink --bfile skatQC_strand_ref --geno 0.01 --hwe 1e-6 --remove plink_ind.txt --make-bed --out ../skatQC_strand_ref_snp --noweb")
system("plink --bfile ../skatQC_strand_ref_snp --recode vcf-iid --allow-no-sex --keep-allele-order --out ../NeuroX_filter_qc")
system("sed -i 's/chr//g' ../NeuroX_filter_qc.vcf")
system("bgzip NeuroX_filter_qc.vcf")
system("tabix -p vcf NeuroX_filter_qc.vcf.gz") 
system("vcftools --gzvcf NeuroX_filter_qc.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out NeuroX_vcftool_qc.flt")


### vcf split
setwd("../")
core <- 9
temp <- fread(input = "NeuroX_vcftool_qc.flt.recode.vcf", sep = "\t", quote = "\n", skip = 30, select = 1:9, header = T, nThread = 9) ### if qc, skip = 3 and select remove

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
fwrite(x = mul_anno_fix, file = "../NeuroX_annotation_fix_0129_rere.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n")

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

## skat-o test
system("mkdir result");setwd("result")
system("pwd")
system("plink --bfile ../skatQC_strand_ref_snp --freq --out skatQC_strand_ref_snp")  ### maf
frq <- fread(file = "skatQC_strand_ref_snp.frq", header = T)
test_fix <- read.table(file = "../NeuroX_annotation_fix_0129_rere.txt", sep = "\t", header = T, stringsAsFactors = F)
# geneset <- read.csv("/home/jinoo/skat-o/parkinson_genset.txt", stringsAsFactors = F,header = F)
geneset <- read.csv("/home/jinoo/skat-o/LSD_geneset.txt", stringsAsFactors = F,header = F)
geneset <- as.character(geneset[,1])

## NULL dataframe
variant_all_parkinson <- data.frame(CHROM = NA, POS = NA, ID = NA, snp131 = NA, REF = NA , ALT = NA, Gene.knownGene = NA, ExonicFunc.knownGene = NA, CADD13_PHRED=NA)[numeric(0), ]

## Parkinson geneset variant
# 1:17312586_G/A   grp format

for(i in 1:length(geneset)){
  variant_all_parkinson <- rbind(variant_all_parkinson, subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")), 
                                                               select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED")))}

# minor allele frequency test
frq$SNP <- gsub(pattern = "chr", replacement = "", x = frq$SNP)
MAF <- mclapply(X = variant_all_parkinson$ID, FUN = function(ID){
  index <- which(x = (str_detect(frq$SNP, paste0("^",ID,"$"))), arr.ind = T)
  return(frq$MAF[index])
}, mc.cores = 9) %>% as.numeric()

variant_all_parkinson <- cbind(variant_all_parkinson, MAF = MAF)

### 1. nonsynonymous geneset
nonsynonymous <- subset(variant_all_parkinson, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" 
                                                | ExonicFunc.knownGene == "stopgain"
                                                | ExonicFunc.knownGene ==  "stoploss" 
                                                | ExonicFunc.knownGene ==  "frameshift_deletion" 
                                                | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                | ExonicFunc.knownGene ==  "frameshift_block_substitution" 
                                                | ExonicFunc.knownGene ==  "splicing"))
nonsynonymous_maf_1 <- subset(nonsynonymous, subset = (MAF <= 0.01), select = "ID")[,1]
nonsynonymous_maf_3 <- subset(nonsynonymous, subset = (MAF <= 0.03), select = "ID")[,1]

### 2. CADD > 12.37 variant
cadd <- subset(variant_all_parkinson, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                 | ExonicFunc.knownGene == "stopgain"
                                                 | ExonicFunc.knownGene ==  "stoploss" 
                                                 | ExonicFunc.knownGene ==  "frameshift_deletion" 
                                                 | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                 | ExonicFunc.knownGene ==  "frameshift_block_substitution" 
                                                 | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
cadd_maf_1 <- subset(cadd, subset = (MAF <= 0.01), select = "ID")[,1]
cadd_maf_3 <- subset(cadd, subset = (MAF <= 0.03), select = "ID")[,1]

### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
lof <- subset(variant_all_parkinson, subset = ( ExonicFunc.knownGene == "stopgain" 
                                                | ExonicFunc.knownGene ==  "stoploss" 
                                                | ExonicFunc.knownGene ==  "frameshift_deletion" 
                                                | ExonicFunc.knownGene ==  "frameshift_insertion" 
                                                | ExonicFunc.knownGene ==  "frameshift_block_substitution" 
                                                | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)

lof_maf_1 <- subset(lof, subset = (MAF <= 0.01), select = "ID")[,1]
lof_maf_3 <- subset(lof, subset = (MAF <= 0.03), select = "ID")[,1]

type_1 <- list(nonsynonymous = nonsynonymous_maf_1, cadd = cadd_maf_1, lof = lof_maf_1)
type_3 <- list(nonsynonymous = nonsynonymous_maf_3, cadd = cadd_maf_3, lof = lof_maf_3)

setID <- list()
for(index in 1:length(type)){
  temp <- data.frame(TYPE=rep(names(type_3[index]), length(type_3[[index]])), stringsAsFactors = F)
  setID[[index]] <- cbind(temp, ID=type_3[[index]])
}
setID <- bind_rows(setID)

system("rm -rf skat.SetID")
write.table(x = setID, "skat_maf_003.SetID", sep = "\t",row.names = F, quote = F, col.names = F)


## skat-o test
Generate_SSD_SetID(File.Bed = "../skatQC_strand_ref_snp.bed",File.Bim = "../skatQC_strand_ref_snp.bim", 
                   File.Fam = "../skatQC_strand_ref_snp.fam", File.SetID = "skat_maf_001.SetID", File.SSD = "skatQC.SSD", File.Info = "skatQC.INFO")
FAM<-Read_Plink_FAM(Filename = "../skatQC_strand_ref_snp.fam", Is.binary = FALSE)
SSD.INFO <- Open_SSD(File.SSD = "skatQC.SSD", File.Info = "skatQC.INFO")

obj<-SKAT_Null_Model(Phenotype ~ Sex, data = FAM, out_type="C", Adjustment = F)
out <- SKAT.SSD.All(SSD.INFO, obj, method = "SKATO")
out$results

### result_adjusted_p value

for(i in 1:length(file)){
  temp <- read.table(file = file[i], header = F)
  system(glue("rm -rf {remove}", remove = file[i]))
  temp_colname <- c("#CHROM","BEGIN","END","MARKER_ID","NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE","STATRHO");colnames(temp) <- temp_colname
  PVALUE_ADJUSTED <- p.adjust(temp$PVALUE, p.adjust.methods[4])
  temp <- cbind(temp, PVALUE_ADJUSTED)
  write.table(x = temp, file = file[i], col.names = T, row.names = F, quote = F, sep = "\t")
}

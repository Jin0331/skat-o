library_load <- function(){
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC);library(tidyverse);library(parallel)
  library(tidyselect);library(magrittr);library(SKAT);
  library(progress);
  library(RMySQL)
  
  # tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/plink/:")
}

# NeuroX missing(REF or ALT unknown )
{
  fread(file = "raw_data/dbGaP_NeuroX_filter.bim", header = F) %>% as_tibble() %>% 
    filter(V5 == '0' | V6 =='0') %>% pull(2) %>% tibble() %>% 
    write_delim(path = "remove_snp_NeuroX.txt", delim = "\t", col_names = F)
  
  ## remove sample 여기에 
  system("plink --bfile raw_data/dbGaP_NeuroX_filter --remove list.id.txt --make-bed --out raw_data/NeuroX_id_rm")
  
  # snp remove
  system("plink --bfile raw_data/NeuroX_id_rm --geno 0.01 --hwe 1e-6 --exclude remove_snp_NeuroX.txt --make-bed --out NeuroX_raw")
  system("./update_build.sh NeuroX_raw NeuroX_15036164_A-b37.Source.strand NeuroX_preQC") ### strand update
  system("plink --bfile NeuroX_preQC --reference-allele NeuroX_15036164_A-b37.strand.RefAlt --make-bed --out NeuroX_preQC_2")
}

# duplicated remove var
# system("plink --bfile NeuroX_preQC_2 --list-duplicate-vars --out dup")
# system("plink --bfile NeuroX_preQC_2 --exclude dup.dupvar --make-bed --out NeuroX_preQC_3")

# sex check =====
if(!dir.exists("2_sexcheck")){
  dir.create("2_sexcheck")
  setwd("2_sexcheck")
} else{
  setwd("2_sexcheck/")
}

# system("plink --bfile ../WES_isec --check-sex --out PPMI")

sex_check <- fread("dbGaP_NeuroX_filter_Gender.sexcheck") %>% as_tibble() %>% 
  select(-IID) %>% mutate(IID = FID) 

check_ <- filter(sex_check, STATUS != "OK")

check_ %>% select(FID, IID) %>% 
  write_delim(path = "remove_sex.txt", delim = "\t")


system("plink --bfile ../NeuroX_preQC_2 --remove remove_sex.txt --make-bed --out ../WES_isec_sex_check")
setwd("../")

# imiss_hetro outlier ====
if(!dir.exists("3_heterozygosity_outliers")){
  dir.create("3_heterozygosity_outliers")
  setwd("3_heterozygosity_outliers")
} else{
  setwd("3_heterozygosity_outliers/")
}

system("plink --bfile ../WES_isec_sex_check --missing --out WES_isec_missing")
system("plink --bfile ../WES_isec_sex_check --het --out WES_isec_het")

plink_het <- fread("WES_isec_het.het", header = T) %>% as_tibble()
plink_het <- plink_het %>% 
  mutate(Heterozygosity_Rate = (`N(NM)` - `E(HOM)`)/`N(NM)`)
plink_imiss <- fread("WES_isec_missing.imiss") %>% as_tibble() %>% select(FID, IID, F_MISS) %>% 
  rename(Proportion_of_missing_genotypes = F_MISS)

pheno <- fread("../NeuroX_preQC_3.fam", header = F) %>% 
  select(V1, V2, V6)
colnames(pheno) <- c("FID","IID","PHENOTYPE")

het_imiss <- left_join(x = plink_het, y = plink_imiss, by = c("FID","IID")) %>% 
  left_join(x = ., y = pheno, by = c("FID", "IID")) %>% 
  mutate(PHENOTYPE = as.factor(PHENOTYPE))

plus_4sd <- mean(het_imiss$Heterozygosity_Rate) + 4 * sd(het_imiss$Heterozygosity_Rate)
minus_4sd <- mean(het_imiss$Heterozygosity_Rate) - 4 * sd(het_imiss$Heterozygosity_Rate)

ggplot(data = het_imiss, aes(x = Proportion_of_missing_genotypes, y = Heterozygosity_Rate, color = PHENOTYPE)) +
  geom_point() +
  geom_vline(xintercept = 0.05, color = "grey") +
  geom_hline(yintercept = minus_4sd, color = "grey")

ggsave("../Hetrozygosity and missing outliers.png")
dev.off()

remove_imiss_hetero <- het_imiss %>%
  filter(Proportion_of_missing_genotypes > 0.05 | (Heterozygosity_Rate < minus_4sd)) %>% 
  select(FID,IID) %>% 
  mutate(IID = FID) %>%
  write_delim(path = "remove_hetero_imiss.txt", delim = "\t")

system("plink --bfile ../WES_isec_sex_check --remove remove_hetero_imiss.txt --make-bed --out ../WES_isec_sex_check_hetero")
setwd("../")

# MDS ====
# Convert vcf to Plink format.
if(!dir.exists("4_MDS")){
  dir.create("4_MDS")
  setwd("4_MDS")
} else{
  setwd("4_MDS/")
}

system("cp ../WES_isec_sex_check_hetero.* .")

### dup var
system("plink --bfile WES_isec_sex_check_hetero --list-duplicate-vars --out dup")
system("plink --bfile WES_isec_sex_check_hetero --exclude dup.dupvar --make-bed --out WES_isec_sex_check_hetero")

NeuroX_bim <- fread("WES_isec_sex_check_hetero.bim") %>% as_tibble()
colnames(NeuroX_bim)[2] <- "ID";colnames(NeuroX_bim)[4] <- "POS"

NeuroX_fix <- fix_load_SKAT(data_name = "NeuroX") %>% 
  select(CHROM, POS, ID, avsnp147)

NeuroX_bim <- left_join(x = NeuroX_bim, y = NeuroX_fix, by = c("ID", "POS"))
NeuroX_bim %>% mutate(ID = ifelse(avsnp147 != ".", avsnp147, ID)) %>% 
  select(-CHROM, -avsnp147) %>% 
  write_delim(path = "WES_isec_sex_check_hetero.bim", delim = "\t", col_names = F)

# 1000G to plink(plink에서 돌릴 것￣)
system("plink --vcf /home/jinoo/skat-o/NeuroX_QC_0804/NeuroX/1000G_1008.vcf.gz --set-missing-var-ids @:#\$1_\$2 --geno 0.2 --allow-no-sex --make-bed --out 1000G_1008")
# Remove individuals based on missing genotype data.
system("plink --bfile 1000G_1008 --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2")
# Remove variants based on MAF.
system("plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3")
# Remove individuals based on missing genotype data.
system("plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4")
# Remove variants based on MAF.
system("plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5")

# NeuroX rs id extraction
temp <- fread(file = "WES_isec_sex_check_hetero.bim", sep = "\t", header = F) %>% 
  as_tibble()

rs_id <- temp$V2 %>% str_extract(pattern = "rs[0-9]+")

for(index in 1:nrow(temp)){
  if(!is.na(rs_id[index]))
    temp$V2[index] <- rs_id[index]
  print(index)
}

write_delim(x = temp, path = "WES_isec_sex_check_hetero.bim", col_names = F)

# Extract the variants present in HapMap dataset from the 1000 genomes dataset.
system("awk '{print$2}' WES_isec_sex_check_hetero.bim > WES_isec_SNPs.txt")
system("plink --bfile 1kG_MDS5 --extract WES_isec_SNPs.txt --make-bed --out 1kG_MDS6")

# Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
system("awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt")
system("plink --bfile WES_isec_sex_check_hetero --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out WES_MDS")


# The datasets now contain the exact same variants.
## The datasets must have the same build. Change the build 1000 Genomes data build., rs157582
system("awk '{print$2,$4}' WES_MDS.map > buildhapmap.txt")

# buildhapmap.txt contains one SNP-id and physical position per line.
system("echo rs157582  > dup.txt")
system("plink --bfile WES_MDS --exclude dup.txt --make-bed --out WES_MDS")
system("plink --bfile 1kG_MDS6 --exclude dup.txt --make-bed --out 1kG_MDS6")

system("plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7")
# 1kG_MDS7 and HapMap_MDS now have the same build.

# 1) set reference genome 
system("awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt")
system("plink --bfile WES_MDS --reference-allele 1kg_ref-list.txt --make-bed --out WES-adj")
# The 1kG_MDS7 and the HapMap-adj have the same reference genome for all SNPs.
# This command will generate some warnings for impossible A1 allele assignment.

# 2) Resolve strand issues.
# Check for potential strand issues.
system("awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp")
system("awk '{print$2,$5,$6}' WES-adj.bim > WES-adj_tmp")
system("sort 1kGMDS7_tmp WES-adj_tmp |uniq -u > all_differences.txt")
# 1624 differences between the files, some of these might be due to strand issues.

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
system("awk '{print$1}' all_differences.txt | sort -u > flip_list.txt")
# Generates a file of 812 SNPs. These are the non-corresponding SNPs between the two files. 
# Flip the 812 non-corresponding SNPs. 
system("plink --bfile WES-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_WES")


# Check for SNPs which are still problematic after they have been flipped.
system("awk '{print$2,$5,$6}' corrected_WES.bim > corrected_WES_tmp")
system("sort 1kGMDS7_tmp corrected_WES_tmp |uniq -u  > uncorresponding_SNPs.txt")
# This file demonstrates that there are 84 differences between the files.

# 3) Remove problematic SNPs from HapMap and 1000 Genomes.
system("awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt")
# The command above generates a list of the 42 SNPs which caused the 84 differences between the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

# Remove the 42 problematic SNPs from both datasets.
system("plink --bfile corrected_WES --exclude SNPs_for_exlusion.txt --make-bed --out WES_MDS2")
system("plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8")

# Merge HapMap with 1000 Genomes Data.
system("plink --bfile WES_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2")

## Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
# Using a set of pruned SNPs
system("plink --bfile MDS_merge2 --genome --out MDS_merge2")
system("plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 20 --out MDS_merge2 --threads 3")

MDS <- fread("MDS_merge2.mds", header = T) %>% as_tibble()
study_fam <- fread("WES-adj.fam") %>% as_tibble() %>% 
  select(FID = V1, IID = V2) %>% 
  mutate(FID = as.character(FID), IID = as.character(IID)) %>% 
  mutate(RACE = "NeuroX")
# system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/20100804.ALL.panel")
G_1000 <- fread("20100804.ALL.panel.txt",sep = "\t") %>% as_tibble() %>% 
  mutate(FID = V1, IID = V1) %>% 
  select(FID, IID, RACE = V2)

# EUR : GBR, FIN, CEU, TSI
# ASN : CHS, CHB, JPT
# AMR : PUR, MXL
# AFR : YRI, LWK, ASW

G_1000$RACE[G_1000$RACE == "GBR"] <- "EUR"
G_1000$RACE[G_1000$RACE == "FIN"] <- "EUR"
G_1000$RACE[G_1000$RACE == "CEU"] <- "EUR"
G_1000$RACE[G_1000$RACE == "TSI"] <- "EUR"

G_1000$RACE[G_1000$RACE == "CHS"] <- "ASN"
G_1000$RACE[G_1000$RACE == "CHB"] <- "ASN"
G_1000$RACE[G_1000$RACE == "JPT"] <- "ASN"

G_1000$RACE[G_1000$RACE == "PUR"] <- "AMR"
G_1000$RACE[G_1000$RACE == "MXL"] <- "AMR"

G_1000$RACE[G_1000$RACE == "YRI"] <- "AFR"
G_1000$RACE[G_1000$RACE == "LWK"] <- "AFR"
G_1000$RACE[G_1000$RACE == "ASW"] <- "AFR"

race_merge <- bind_rows(study_fam, G_1000) %>% 
  mutate(RACE = as.factor(RACE))

MDS <- MDS %>% left_join(x = ., y = race_merge, by = c("FID", "IID"))

# MDS <- fread("WES_QC_0821.mds", header = T)
### plot
ggplot(data = MDS, aes(x = C1, y = C2, color = RACE)) +
  geom_point() +
  ggtitle("NeuroX") + 
  labs(x = "1 component", y = "2 component") + 
  theme_classic() + 
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 25, color = "black"),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(face = "bold", size = 10),
        legend.box.background = element_rect(fill = "black", size = 1.1),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 12))
ggsave("../NeuroX_mds_1_2.png")
dev.off()

### QC
WES_MDS <- MDS %>% filter(RACE == "NeuroX") %>% 
  select(-RACE)
col_mean <- apply(WES_MDS[,4:23], 2, mean)
col_sd <- apply(WES_MDS[,4:23], 2, sd)

removal_MDS <- c()
for(i in 1:nrow(WES_MDS)){
  for( j in 4:7){ 
    if( (WES_MDS[i,j] > (col_mean[j-3] + 4 * col_sd[j-3])) | (WES_MDS[i,j] < (col_mean[j-3] - 4 * col_sd[j-3]))){
      removal_MDS <- c(removal_MDS, WES_MDS[i,1])}
  }
}

removal_MDS %>% as.character() %>% unique() %>% as_tibble() %>% 
  mutate(FID = value, IID = value) %>% select(-value) %>% 
  write_delim(x = ., path = "remove_mds.txt", delim = "\t")

#### mds result
# Extract these individuals in HapMap data.
system("plink --bfile ../WES_isec_sex_check_hetero --remove remove_mds.txt --make-bed --out ../WES_QC_0821")
system("plink --bfile ../WES_QC_0821 --extract WES_isec_sex_check_hetero.prune.in --genome --out ../WES_QC_0821")
system("plink --bfile ../WES_QC_0821 --read-genome ../WES_QC_0821.genome --cluster --mds-plot 20 --out ../WES_QC_0821")


setwd("../")
system("plink --bfile WES_QC_0821 --indep-pairwise 50 5 0.2 --out skat_pruning")
system("plink --bfile WES_QC_0821 --genome full --out skatQC_IBD") ### MDS
system("plink --bfile WES_QC_0821 --read-genome skatQC_IBD.genome --extract skat_pruning.prune.in --mds-plot 20 --cluster --out NeuroX_0821")


system("plink --bfile WES_QC_0821 --freq --out NeuroX_0821")

# making cov

mds <- fread(file = "NeuroX_0821.mds", header = T) %>% 
  select(-IID) %>% 
  mutate(IID = FID)

age <- fread(file = "NeuroX_phenotype.txt", sep = "\t") %>% as_tibble() %>% 
  select(IID = SUBJECT_ID, AGE)

Fmiss <- fread(file = "NeuroX_0821.imiss") %>% as_tibble() %>% 
  select(IID, F_MISS)


cov <- left_join(x = mds, y = age, by = c("IID")) %>% 
  left_join(x = ., y = Fmiss, by = c("IID")) %>% 
  as_tibble() %>% 
  select(FID, IID, AGE, F_MISS,	C1,	C2,	C3,	C4)

fwrite(x = cov, file ="NeuroX_0821.cov", row.names = F, sep = "\t")



# removal_MDS %>% as.character() %>% unique() %>% as_tibble() %>% 
#   mutate(FID = value, IID = value) %>% select(-value) %>% 
#   write_delim(x = ., path = "remove_mds.txt", delim = "\t")
# WES_MDS %>% select(-IID, -SOL) %>% 
#   write_delim(x = ., path = "WES_MDS.txt", delim = "\t")

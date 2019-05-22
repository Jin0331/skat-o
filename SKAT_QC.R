# library path
{
  library(glue);library(vcfR);library(data.table);library(foreach);library(doMC)
  library(tidyverse)
# tool path 
  Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/jinoo/tool/:")
}

## unknown snp QC(e.g. minor | major 0) for strand
{
  system("mkdir data");setwd("data/")
  
  file <- list.files(path = "../", pattern = ".bim", full.names = T)bim_qc <- fread(file = file, header = F, stringsAsFactors = F, data.table = F)
  bim_zero_snp <- subset.data.frame(x = bim_qc, subset = (V5 == '0' | V6 =='0'), select = "V2")[,1]
  write.table(x = bim_zero_snp, file = "remove_snp_file.txt", sep = "\t", quote = F, row.names = F, col.names = F)
  system("plink --bfile ../pd_v3_550v3 --exclude remove_snp_file.txt --make-bed --out skatQC")
  system("./update_build.sh skatQC HumanHap550v3_A-b36.Ilmn.strand skatQC_strand") ### strand update
  system("plink --bfile skatQC_strand --reference-allele HumanHap550v3_A.b36.RefAlt --chr 1-22, 23-24 --make-bed --out skatQC_strand_ref")
  setwd("../")
}

## sample QC
{
  #QC
  {
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
          removal_MDS <- c(removal_MDS, plink_MDS[i,2])
      }
    }
    
    #3 missing
    plink_missing <- read.table(file = "skatQC_missing.imiss", header = T, stringsAsFactors = F)
    removal_imissing <- subset(plink_missing, subset = (F_MISS > 0.15), select = c("IID","F_MISS"))[,1]
    
    #4 gender
    plink_gender <- read.table(file = "skatQC_sex.sexcheck", header = T, stringsAsFactors = F)
    removal_gender <- subset(plink_gender, subset = (STATUS == "PROBLEM"), select = c("IID"))[,1]
    
    #5 ajkvalue- skip
    {
    #5 vcf QC -- vcftools, ajk value
    
    # system("pseq NeuroX new-project --resources /home/lee/pseq/hg19/")
    # system("pseq NeuroX.pseq load-plink --file skatQC_snp --id skatQC_snp")
    # system("pseq NeuroX.pseq write-vcf > NeuroX_filter.vcf")
    # system("pseq NeuroX_filter.vcf write-vcf --format BGZF --file NeuroX_filter.vcf.gz")
    # system("mv NeuroX_filter.vcf NeuroX_filter.vcf.bak")
    # # system("gzip -d NeuroX_filter.vcf.gz")
    # 
    # system("plink --bfile ../data/skatQC_strand_ref --recode vcf-iid --allow-no-sex --keep-allele-order --out IPDGC_filter")
    # system("sed -i 's/chr//g' IPDGC_filter.vcf")
    # 
    }
    
    #6 merge & write
    removal_iid <- unique(c(removal_MDS,removal_het,removal_imissing, removal_gender)) 
    write.table(x = removal_iid, file = "plink_ind.txt",sep = "\n", quote = F, row.names = F, col.names = F)
    rm(list=ls());gc()
    
    system("plink --bfile ../data/skatQC_strand_ref --geno 0.01 --hwe 1e-6 --make-bed --out skatQC_strand_ref_snp --noweb")
    system("plink --bfile skatQC_strand_ref_snp --recode vcf-iid --allow-no-sex --keep-allele-order --out ../IPDGC_GWAS_filter_qc")
    system("sed -i 's/chr//g' ../IPDGC_GWAS_filter_qc.vcf")
    system("bgzip ../IPDGC_GWAS_filter_qc.vcf")
    system("tabix -p vcf ../IPDGC_GWAS_filter_qc.vcf.gz") 
    
    setwd("../")
    system("vcftools --gzvcf IPDGC_GWAS_filter_qc.vcf.gz --min-alleles 2 --max-alleles 2 --remove plink/plink_ind.txt --recode --out IPDGC_qc_com")
  }
  
  #Annotation
  {
    ### vcf split
    system("mkdir anno");setwd("anno/")
    core <- 9
    temp_tidy <- fread(input = "../IPDGC_qc_com.recode.vcf", sep = "\t", quote = "\n", skip = 29, header = T, nThread = 5) %>% as_tibble()
    
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
    registerDoMC(9)
    system.time(
      foreach(i=1:core) %dopar% {
        system(glue("table_annovar.pl {path[i]} /home/lee/annovar/humandb/ -buildver hg18 -out {path[i]} -remove -protocol knownGene,dbnsfp30a -operation g,f -nastring . -vcfinput"))
      }
    );gc() #
    
    #### make annotated fix
    system("rm -rf *.txt");system("rm -rf *.avinput")
    temp <- list.files(pattern = "hg18")
    
    mul_vcf <- mclapply(temp, read.vcfR, convertNA =T, checkFile = F, mc.cores = 2)
    mul_anno_fix <- lapply(X = mul_vcf, FUN = function(temp){
      temp2 <- vcfR2tidy(temp, info_only = T, single_frame = F, toss_INFO_column = T)
      return(temp2$fix)
    })
    mul_anno_fix <- bind_rows(mul_anno_fix)
    rm(mul_vcf)
    fwrite(x = mul_anno_fix, file = "../IPDGC_GWAS_fix.txt", quote = F, sep = "\t", row.names = F, col.names = T, eol = "\n")
    
    setwd("../")

  }
  
  ## new ped
  setwd("../")
  skat_row_ped <- read.table(file = "../NeuroX_filter_epact_row.ped", header = T, stringsAsFactors = F)
  skat_subset <- read.table(file = "plink_ind.txt", header = F, stringsAsFactors = F)[,1]
  '%!in%' <- function(x,y)!('%in%'(x,y))
  skat_QC_ped <- NULL
  skat_QC_ped <- subset(skat_row_ped, subset = skat_row_ped$fid %!in% skat_subset)
  write.table(x = skat_QC_ped, file = "../NeuroX_epacts_0129.ped", col.names = T, row.names = F, sep = "\t", quote = F)

}

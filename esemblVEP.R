BiocManager::install("ensemblVEP", version = "3.8")
library(ensemblVEP)

system("zcat /home/jinoo/skat-o/SKAT_data/IPDGC_QC.vcf.gz | cut -f1-9 > IPDGC.vars.vcf")
system("bgzip IPDGC.vars.vcf")
temp_vcf <-file.path("/home/jinoo/skat-o/vep_anno/IPDGC.vars.vcf.gz")

param <- VEPFlags(flags = list(sift="b", polyphen = "b", vcf = TRUE, 
                      cache = TRUE, dir = "/home/jinoo/ensembl-vep/data/",
                      format = "vcf", assembly = "GRCh37",output_file = "/home/jinoo/skat-o/vep_anno/IPDGC_vep.vcf", offline = TRUE,
                      force_overwrite = TRUE, pick = TRUE, verbose = TRUE))

param
flags(param)
ensemblVEP(temp_vcf, param)



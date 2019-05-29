library(dplyr)
library(stringr)

Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/lee/tool/:/home/lee/annovar/:/home/lee/epacts_0913/bin/")

#### geneset
geneset_merge <- read.csv(file = "../1 gene sets 190208 for SKAT_O.csv", header = T, encoding = "UTF-8", stringsAsFactors = F)
col_name <- c("GENE","O1","O2","O3","O4","O5","O6","O7","O8","O9", "M1", "M2")
colnames(geneset_merge) <- col_name

M1 <- subset.data.frame(x = geneset_merge, subset = (str_detect(string = M1, pattern = "[^0]")), select = "GENE")[,1]
M1 <- c(M1, rep(x = "", (nrow(geneset_merge) - length(M1))))

M2 <- subset.data.frame(x = geneset_merge, subset = (str_detect(string = M2, pattern = "[^0]")), select = "GENE")[,1]
M2 <- c(M2, rep(x = "", (nrow(geneset_merge) - length(M2))))

geneset_merge$M1 <- M1;geneset_merge$M2 <- M2

########## M3, M4 ###########
C <- unique.default(geneset_merge$O2);D <- unique.default(geneset_merge$O1) ## M3
A <- unique.default(geneset_merge$M1);B <- unique.default(geneset_merge$O1) ## M4
M3 <- C[!(C%in%D)];M3 <- c(M3, rep(x = "", (nrow(geneset_merge) - length(M3))))
M4 <- A[!(A%in%B)];M4 <- c(M4, rep(x = "", (nrow(geneset_merge) - length(M4))))
geneset_merge <- cbind(geneset_merge, M3);geneset_merge <- cbind(geneset_merge, M4)
col_name <- colnames(geneset_merge)
geneset_merge$M3 <- as.character(geneset_merge$M3);geneset_merge$M4 <- as.character(geneset_merge$M4)
rm(A,B,C,col_name,D,M1,M2,M3,M4)
#############################


#### non-syn, cadd, lof
system("mkdir result");setwd("result");system("pwd")
path <- list.files(path = "/home/jinoo/skat-o/0208_dataset/", pattern = "(fix)", full.names = T)
dataset_name <- c("IPDGC","NeuroX","PPMI");col_name <- colnames(geneset_merge)

for(index in 1:1){ #1:length(dataset_name
  print(dataset_name)
  test_fix <- read.table(file = path[2], sep = "\t", header = T, stringsAsFactors = F)
  
  for(colname in c(2)){ #for(colname in 1:ncol(geneset_merge))
    geneset <- grep(x = geneset_merge[,colname], pattern = "[^ ]", value = T)
    print(col_name[colname])
    variant_all_parkinson <- data.frame(CHROM = NA, POS = NA, ID = NA, snp131 = NA, REF = NA , ALT = NA, Gene.knownGene = NA, ExonicFunc.knownGene = NA, CADD13_PHRED=NA)[numeric(0), ]
    for(i in 1:length(geneset)){
      variant_all_parkinson <- rbind(variant_all_parkinson, subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")), 
                                                                   select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED")))
    }
    
    ### 1. nonsynonymous geneset
    nonsynonymous <- subset(variant_all_parkinson, subset = ((ExonicFunc.knownGene == "nonsynonymous_SNV")) | ExonicFunc.knownGene == "stopgain" 
                            | ExonicFunc.knownGene ==  "stoploss" | ExonicFunc.knownGene ==  "frameshift_deletion"
                            | ExonicFunc.knownGene ==  "frameshift_insertion" | ExonicFunc.knownGene ==  "frameshift_block_substitution" 
                            | ExonicFunc.knownGene ==  "splicing")
    
    i <- NULL;nonsynonymous_geneset <- c()
    for(i in 1:nrow(nonsynonymous)){
      nonsynonymous_geneset <- c(nonsynonymous_geneset, paste0(nonsynonymous[i,1], ":", nonsynonymous[i,2], "_", nonsynonymous[i,4], "/", nonsynonymous[i,5]))
      print(i)
    }
    
    ### 2. CADD > 12.37 variant
    cadd <- subset(variant_all_parkinson, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                     | ExonicFunc.knownGene == "stopgain"
                                                     | ExonicFunc.knownGene ==  "stoploss" 
                                                     | ExonicFunc.knownGene ==  "frameshift_deletion" 
                                                     | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                     | ExonicFunc.knownGene ==  "frameshift_block_substitution" 
                                                     | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
    
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
    nonsynonymous_geneset <- t(data.frame(nonsynonymous_geneset, stringsAsFactors = F));rownames(nonsynonymous_geneset) <- paste0("nonsynonymous","_",col_name[colname])
    cadd_geneset <- t(data.frame(cadd_geneset, stringsAsFactors = F));rownames(cadd_geneset) <- paste0("cadd","_",col_name[colname])
    if(nrow(lof) !=0){
      lof_geneset <- t(data.frame(lof_geneset, stringsAsFactors = F));rownames(lof_geneset) <- paste0("lof","_",col_name[colname])
    }
    print(length(nonsynonymous_geneset))
    print(length(cadd_geneset))
    print(length(lof_geneset))
    #system("rm -rf skat_grp.grp")
    write.table(x = nonsynonymous_geneset, file = paste0(dataset_name[index], ".grp"), sep = "\t",row.names = T, quote = F, col.names = F, append = T)
    write.table(x = cadd_geneset, file = paste0(dataset_name[index], ".grp"), sep = "\t", row.names = T, quote = F, col.names = F, append = T)
    write.table(x = lof_geneset, file = paste0(dataset_name[index], ".grp"), sep = "\t", row.names = T, quote = F, col.names = F, append = T)
    
    # rm(cadd, cadd_geneset, lof, lof_geneset, nonsynonymous, nonsynonymous_geneset, geneset, variant_all_parkinson)
  }
  
   ### gene
  # variant_gene <- unique(geneset_merge$GENE)
  # paste_temp <- list()
  # temp <- NULL
  # for( i in 1:length(variant_gene)){
  #   temp <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED > 12.37))
  #   temp <- subset.data.frame(temp, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"| ExonicFunc.knownGene == "stopgain"
  #                  | ExonicFunc.knownGene ==  "stoploss" | ExonicFunc.knownGene ==  "frameshift_deletion"
  #                  | ExonicFunc.knownGene ==  "frameshift_insertion" | ExonicFunc.knownGene ==  "frameshift_block_substitution"
  #                  | ExonicFunc.knownGene ==  "splicing"))
  #   if(nrow(temp) != 0){
  #     for(j in 1:nrow(temp)){
  #       paste_temp[[j]] <- paste0(temp[j,1], ":", temp[j,2],"_", temp[j,4],"/",temp[j,5])
  #     }
  #     paste_temp <- t(data.frame(as.character(paste_temp),stringsAsFactors = F));rownames(paste_temp) <- variant_gene[i]
  #     write.table(x = paste_temp, file = paste0(dataset_name[index], "_gene.grp"), sep = "\t", row.names = T, col.names = F, append = T, quote = F)
  #   }
  #   print(i)
  #   paste_temp <- list()
  # }
  # rm(test_fix, paste_temp, temp)
}# 1:17312586_G/A   grp format



### EPACTS skat
MAF <- c(0.01, 0.03)
for(i in MAF){
  # print(i);print("IPDGC")
  # system(glue("epacts-group -vcf /home/jinoo/skat-o/0208_dataset/IPDGC_QC.vcf.gz -groupf IPDGC.grp  -out IPDGC_maf{maf}.skat -ped /home/jinoo/skat-o/0208_dataset/IPDGC.ped -max-maf {maf} -pheno disease -cov sex -missing ./. -test skat -skat-o -run 8", maf = i))
  # print("PPMI")
  # system(glue("epacts-group -vcf ../PPMI_QC.vcf.gz -groupf PPMI.grp  -out PPMI_maf{maf}.skat -ped ../PPMI.ped -max-maf {maf} -pheno disease -cov sex -missing ./. -test skat -skat-o -run 8", maf = i))
  # print("NeuroX")
  system(glue("epacts-group -vcf /home/jinoo/skat-o/0208_dataset/NeuroX_QC.vcf.gz -groupf NeuroX.grp  -out NeuroX_maf{maf}.skat -ped /home/jinoo/skat-o/0208_dataset/NeuroX.ped -max-maf {maf} -pheno disease -cov sex -missing ./. -test skat -skat-o -run 2", maf = i))
}

### result_adjusted_p value
# system("find . ! -name '*.epacts' -delete")
file <- list.files(pattern = "(.epacts)")

for(i in 1:length(file)){
  temp <- read.table(file = file[i], header = F)
  system(glue("rm -rf {remove}", remove = file[i]))
  temp_colname <- c("#CHROM","BEGIN","END","MARKER_ID","NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE","STATRHO");colnames(temp) <- temp_colname
  PVALUE_ADJUSTED <- p.adjust(temp$PVALUE, p.adjust.methods[4])
  temp <- cbind(temp, PVALUE_ADJUSTED)
  write.table(x = temp, file = file[i], col.names = T, row.names = F, quote = F, sep = "\t")
}



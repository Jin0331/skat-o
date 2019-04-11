library(SKAT);library(tidyverse);library(data.table)
library(parallel);library(doParallel)


# function
geneset_load <- function(){
  return_list <- list()
  geneset_merge <- fread(file = "/home/jinoo/skat-o/SKAT_data/1 gene sets 190308 for SKAT_O.txt", header = T, 
                         sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
  
  col_name <- c("GENE","O1","O2","O2M1","O2M2","O2M3","O3","O4","O5","O6", "O7", "O8", "O9","M2","M2P")
  colnames(geneset_merge) <- col_name
  
  one_hot <- lapply(X=select(geneset_merge,-GENE, -M2), FUN = function(x){ifelse(nchar(x) > 0, 1, 0)}) %>% bind_cols()
  geneset_onehot <- cbind(GENE=geneset_merge$GENE, one_hot) %>% mutate(GENE = as.character(GENE)) %>% cbind(M2 = geneset_merge$M2)
  
  for(geneset in 2:ncol(geneset_onehot)){
    for(k in 1:nrow(geneset_onehot)){
      if(geneset_onehot[k, geneset] >= 1){
        geneset_onehot[k,geneset] = geneset_onehot[k,1]
      } else{
        geneset_onehot[k,geneset] = NA
      }
      
    }
  }
  return_list[[1]] <- geneset_onehot %>% as.list()
  return_list[[2]]<- c("GENE","O1","O2","O2M1","O2M2","O2M3","O3","O4","O5","O6", "O7", "O8", "O9","M2","M2P", "EPI","CHD")
  
  ## negative geneset
  temp <- fread(file = "/home/jinoo/skat-o/SKAT_data/denovo_chd_epi.txt", header = T, sep = "\t")
  return_list[[1]][[16]] <-  temp[[1]]
  return_list[[1]][[17]] <-  temp[[2]][temp[[2]] != ""]
  
  return(return_list)
} ## 1 = geneset_list, 2 = col_name

fix_load <- function(data_name){
  test_fix <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name, "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
    as_tibble()
  test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
  test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
  test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
  
  test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
  for(change in 1:length(test_id)){
    if(test_id[change] == ""){
      test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
    }
  }
  test_fix$ID <- test_id
  
  return(test_fix)
} # data_name, "IPDGC", "NeuroX"

geneset_setid <- function(geneset_merge, col_name, test_fix, data_name, index_){
  for(gene in index_){ # 2:length(col_name), c(2,3,7,16,17)
    geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
    print(col_name[gene])
    variant <- list()
    ## Parkinson geneset variant
    for(i in 1:length(geneset)){
      variant[[i]] <- subset(test_fix, subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")),
                             select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","CADD13_PHRED"))
    }
    variant <- bind_rows(variant)
    
    ### 1. nonsynonymous geneset
    nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                 | ExonicFunc.knownGene == "stopgain"
                                                 | ExonicFunc.knownGene ==  "stoploss"
                                                 | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                 | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                 | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                 | ExonicFunc.knownGene ==  "splicing"
      ))
    nonsynonymous <- subset(nonsynonymous, select = "ID")[,1]
    
    
    ### 2. CADD > 12.37 variant
    cadd <- subset(variant, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                         | ExonicFunc.knownGene == "stopgain"
                                         | ExonicFunc.knownGene ==  "stoploss"
                                         | ExonicFunc.knownGene ==  "frameshift_deletion"
                                         | ExonicFunc.knownGene ==  "frameshift_insertion"
                                         | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                         | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
    cadd <- subset(cadd, select = "ID")[,1]
    
    ### 3. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
    lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                        | ExonicFunc.knownGene ==  "stoploss"
                                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                        | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
      
    lof <- subset(lof, select = "ID")[,1]

    type <- list(nonsynonymous = nonsynonymous, cadd = cadd, lof = lof)
      
    setID <- list()
    for(j in 1:length(type)){
      temp <- data.frame(TYPE=rep(names(type[j]), length(type[[j]])), stringsAsFactors = F)
      setID[[j]] <- cbind(temp, ID=type[[j]])
    }
    setID <- bind_rows(setID)
    setID$TYPE <- paste0(setID$TYPE, "__",col_name[gene])
    # system("rm -rf skat.SetID")
    fwrite(x = setID, file = paste0(data_name, "_skat.SetID"), sep = "\t",row.names = F, quote = F, col.names = F,append = T)
  }
} 

gene_setid <- function(all_gene, test_fix, data_name){ 
  variant_gene <- unique(all_gene)
  temp <- list()
  for(i in 1:length(variant_gene)){
    t1 <- subset.data.frame(test_fix, subset = (Gene.knownGene %in% variant_gene[i] & CADD13_PHRED > 12.37))
    temp[[i]] <- subset.data.frame(t1, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV" ), 
                                   select = c("Gene.knownGene","ID","ID2", "CADD13_PHRED"))
  }
  setID <- bind_rows(temp)
  
  write.table(x = select(setID, Gene.knownGene, ID), file = paste0(data_name, "_skat_gene.SetID"), sep = "\t",row.names = F, quote = F, col.names = F,append = T)
  write.table(x = select(setID, Gene.knownGene, ID2, CADD13_PHRED), file = paste0(data_name, "_skat_gene_variant.txt"), 
              sep = "\t",row.names = F, quote = F, col.names = F,append = T)
} # all_gene, column 1 of geneset

{
#   run_skat_gender_age <- function(data_name, flag = "default"){
#   
#   if(flag == "geneset"){
# 
#       Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
#                          File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
#                          File.Fam = paste0("../SKAT_data/",data_name,".fam"),
#                          File.SetID = paste0(data_name,"_skat.SetID"), 
#                          File.SSD = paste0(data_name,".SSD"), 
#                          File.Info = paste0(data_name,".INFO"))
#       
#       
#       FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
#                                 File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
#       SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
#       
#       # if(length(colnames(FAM)) > 20){
#       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + 
#       #                           C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
#       #                           C19 + C20, data = FAM, out_type="C", Adjustment = F)
#       # }else{
#       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4, 
#       #                         data = FAM, out_type="C", Adjustment = F)
#       # }
#       
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE, data = FAM, out_type="C", Adjustment = F)
#     
#       ## run
#       out <- list()
#       out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
#       out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
#         
#       fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_result.txt"), col.names = T, row.names = F, sep = "\t")
#       fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_result.txt"), col.names = T, row.names = F, sep = "\t")
#       Close_SSD()
#     
#   } else {
#       Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
#                          File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
#                          File.Fam = paste0("../SKAT_data/",data_name,".fam"),
#                          File.SetID = paste0(data_name,"_skat_gene.SetID"), 
#                          File.SSD = paste0(data_name,"_gene.SSD"), 
#                          File.Info = paste0(data_name,"_gene.INFO"))
#       
#       
#       FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
#                                 File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
#       SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
#       
#       # if(length(colnames(FAM)) > 20){
#       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + 
#       #                           C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
#       #                           C19 + C20, data = FAM, out_type="C", Adjustment = F)
#       # }else{
#       #   obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4, 
#       #                         data = FAM, out_type="C", Adjustment = F)
#       # }
#       
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE, data = FAM, out_type="C", Adjustment = F)
#       
#       ## run
#       out <- list()
#       out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
#       out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
#       
#       fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result.txt"), col.names = T, row.names = F, sep = "\t")
#       fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result.txt"), col.names = T, row.names = F, sep = "\t")
#       Close_SSD()
#       
#     }
# 
# } # flag (geneset, gene)
#   run_skat_all_SSD_ALL <- function(data_name, flag = "default"){
#   
#   if(flag == "geneset"){
#     
#     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
#                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
#                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
#                        File.SetID = paste0(data_name,"_skat.SetID"), 
#                        File.SSD = paste0(data_name,".SSD"), 
#                        File.Info = paste0(data_name,".INFO"))
#     
#     
#     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
#                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
#     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
#     
#     if(length(colnames(FAM)) > 20){
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
#                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
#                               C19 + C20, data = FAM, out_type="C", Adjustment = F)
#     }else{
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
#                             data = FAM, out_type="C", Adjustment = F)
#     }
#     
#     
#     ## run
#     out <- list()
#     out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
#     out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
#     
#     fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
#     fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
#     Close_SSD()
#     
#   } else {
#     Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
#                        File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
#                        File.Fam = paste0("../SKAT_data/",data_name,".fam"),
#                        File.SetID = paste0(data_name,"_skat_gene.SetID"), 
#                        File.SSD = paste0(data_name,"_gene.SSD"), 
#                        File.Info = paste0(data_name,"_gene.INFO"))
#     
#     
#     FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
#                               File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
#     SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
#     
#     if(length(colnames(FAM)) > 20){
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
#                               C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
#                               C19 + C20, data = FAM, out_type="C", Adjustment = F)
#     }else{
#       obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
#                             data = FAM, out_type="C", Adjustment = F)
#     }
#     
#     
#     ## run
#     out <- list()
#     out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
#     out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
#     
#     fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
#     fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
#     Close_SSD()
#   }
#   
# }
}

run_skat_all_cov <- function(data_name, flag = "default"){
  
  if(flag == "geneset"){
    
    Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
                       File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
                       File.Fam = paste0("../SKAT_data/",data_name,".fam"),
                       File.SetID = paste0(data_name,"_skat.SetID"), 
                       File.SSD = paste0(data_name,".SSD"), 
                       File.Info = paste0(data_name,".INFO"))
    
    
    FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                              File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))

    # skat run
    {
      cl <- makeCluster(detectCores() - 1)
      # clusterExport(cl, "obj")
      clusterExport(cl, "data_name")
      # clusterExport(cl, "SSD.INFO")
      # clusterEvalQ(cl, obj)
      clusterEvalQ(cl, data_name)
      clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
    }
    
    SKAT_result <- 
      parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index){
        FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                                  File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE)
        SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
        
        if(data_name == "IPDGC"){
          obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                                  C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                                  C19 + C20 + F_MISS, data = FAM, out_type="C", Adjustment = T)
        }else{
          obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + F_MISS,
                                data = FAM, out_type="C", Adjustment = F)
        }
        
        Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
        out_03 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03)
        out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01)

        skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker,
               n_marker_test = out_03$param$n.marker.test, p_value = out_03$p.value)
        skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker,
               n_marker_test = out_01$param$n.marker.test, p_value = out_01$p.value)

        return(bind_cols(skat_001,skat_003))
        
      })
    
    stopCluster(cl)
    
    results <- SKAT_result %>% bind_rows()
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all_NR.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
    
    
    
  } else {
    Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
                       File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
                       File.Fam = paste0("../SKAT_data/",data_name,".fam"),
                       File.SetID = paste0(data_name,"_skat_gene.SetID"), 
                       File.SSD = paste0(data_name,"_gene.SSD"), 
                       File.Info = paste0(data_name,"_gene.INFO"))
    
    
    FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                              File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
    
    # skat run
    {
      cl <- makeCluster(detectCores() - 1)
      # clusterExport(cl, "obj")
      clusterExport(cl, "data_name")
      # clusterExport(cl, "SSD.INFO")
      # clusterEvalQ(cl, obj)
      clusterEvalQ(cl, data_name)
      clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
    }
    
    SKAT_result <- 
      parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index){
        FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                                  File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE)
        SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
        
        if(data_name == "IPDGC"){
          obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                                  C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                                  C19 + C20 + F_MISS, data = FAM, out_type="C", Adjustment = T)
        }else{
          obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + F_MISS,
                                data = FAM, out_type="C", Adjustment = F)
        }
        
        Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
        out_03 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03)
        out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01)
        
        skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker,
                           n_marker_test = out_03$param$n.marker.test, p_value = out_03$p.value)
        skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker,
                           n_marker_test = out_01$param$n.marker.test, p_value = out_01$p.value)
        
        return(bind_cols(skat_001,skat_003))
        
      })
    
    stopCluster(cl)
    
    results <- SKAT_result %>% bind_rows()
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all_NR_gene.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
  }
  
}
run_skat_all_resampling <- function(data_name, flag = "default"){
  
  if(flag == "geneset"){
    
    Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
                       File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
                       File.Fam = paste0("../SKAT_data/",data_name,".fam"),
                       File.SetID = paste0(data_name,"_skat.SetID"), 
                       File.SSD = paste0(data_name,".SSD"), 
                       File.Info = paste0(data_name,".INFO"))
    
    
    FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                              File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,".SSD"), File.Info = paste0(data_name,".INFO"))
    
    if(length(colnames(FAM)) > 20){
      obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                              C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                              C19 + C20, data = FAM, out_type="C", Adjustment = T, n.Resampling = 5000)
    }else{
      obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
                            data = FAM, out_type="C", Adjustment = F, n.Resampling = 5000)
    }
    
    # skat run
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, "obj")
    # clusterExport(cl, "SSD.INFO")
    # clusterEvalQ(cl, obj)
    clusterEvalQ(cl, {library(SKAT);library(tidyverse);library(data.table)})
    
    
    SKAT_result <- 
      parLapply(cl = cl, X = 1:nrow(SSD.INFO$SetInfo), fun = function(SET_index, name = "IPDGC"){
        
        SSD.INFO <- Open_SSD(File.SSD = paste0(name,".SSD"), File.Info = paste0(name,".INFO"))
        Z <- Get_Genotypes_SSD(SSD_INFO = SSD.INFO, SET_index, is_ID = T)
        out_03 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.03)
        out_01 <- SKAT(Z, obj, method = "optimal.adj", missing_cutoff = 0.9, max_maf = 0.01)
        resample_p_03 <- Get_Resampling_Pvalue(out_03)[[1]]
        resample_p_01 <- Get_Resampling_Pvalue(out_01)[[1]]
        
        skat_003 <- tibble(SetID_003 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_03$param$n.marker, 
                           n_marker_test = out_03$param$n.marker.test, 
                           emprical_pvalue = resample_p_03, p_value = out_03$p.value)
        skat_001 <- tibble(SetID_001 = SSD.INFO$SetInfo[SET_index,2], n_marker = out_01$param$n.marker, 
                           n_marker_test = out_01$param$n.marker.test, 
                           emprical_pvalue = resample_p_01, p_value = out_01$p.value)
        
        return(bind_cols(skat_001,skat_003))
        
      })
    
    stopCluster(cl)
    
    results <- SKAT_result %>% bind_rows()
    
    fwrite(x = results, file = paste0(data_name,"_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
    
    
    
  } else {
    Generate_SSD_SetID(File.Bed = paste0("../SKAT_data/",data_name,".bed"),
                       File.Bim = paste0("../SKAT_data/",data_name,".bim"), 
                       File.Fam = paste0("../SKAT_data/",data_name,".fam"),
                       File.SetID = paste0(data_name,"_skat_gene.SetID"), 
                       File.SSD = paste0(data_name,"_gene.SSD"), 
                       File.Info = paste0(data_name,"_gene.INFO"))
    
    
    FAM <- Read_Plink_FAM_Cov(Filename = paste0("../SKAT_data/",data_name,".fam"),
                              File_Cov = paste0("../SKAT_data/",data_name,".cov"), Is.binary = FALSE )
    SSD.INFO <- Open_SSD(File.SSD = paste0(data_name,"_gene.SSD"), File.Info = paste0(data_name,"_gene.INFO"))
    
    if(length(colnames(FAM)) > 20){
      obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 +
                              C9 + C10 + C11 + C12 + C13 + C14 + C15 + C16 + C17 + C18 +
                              C19 + C20, data = FAM, out_type="C", Adjustment = F)
    }else{
      obj <-SKAT_Null_Model(Phenotype ~ Sex + AGE + C1 + C2 + C3 + C4,
                            data = FAM, out_type="C", Adjustment = F)
    }
    
    
    ## run
    out <- list()
    out[[1]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.01, missing_cutoff = 0.9)
    out[[2]] <- SKAT.SSD.All(SSD.INFO = SSD.INFO, obj, method = "optimal.adj", max_maf = 0.03, missing_cutoff = 0.9)
    
    fwrite(x = out[[1]]$results, file = paste0(data_name,"_0.01_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
    fwrite(x = out[[2]]$results, file = paste0(data_name,"_0.03_gene_result_cov_all.txt"), col.names = T, row.names = F, sep = "\t")
    Close_SSD()
  }
  
}


### cov making
{
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --indep-pairwise 50 5 0.2 --out skat_pruning", name = "NeuroX"))
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --genome full --out skatQC_IBD", name = "NeuroX")) 
  # system(glue("plink --bfile /home/jinoo/skat-o/SKAT_data/{name} --read-genome skatQC_IBD.genome --extract skat_pruning.prune.in --mds-plot 20 --cluster --out {name}_MDS",
  #             name = "NeuroX"))
  # 
  # 
  # ###
  # 
  # FAM <- Read_Plink_FAM_Cov(Filename = "../IPDGC.fam", File_Cov = "IPDGC.cov",Is.binary = F)
  # 
  # IPDGC_mds <- fread(file = "IPDGC_MDS.mds", header = T) %>% select(., -SOL)
  # IPDGC_age_control <- fread(file = "IPDGC_control_phenotype.txt", select = c(2,7)) %>% select(., IID = SUBJECT_ID, AGE = Age)
  # IPDGC_age_case <- fread(file = "IPDGC_case_phenotype.txt", select = c(2,7)) %>% select(., IID = SUBJECT_ID, AGE = Age)
  # IPDGC_age <- rbind(IPDGC_age_case, IPDGC_age_control)
  # 
  # IPDGC_cov <- left_join(x = FAM, y = IPDGC_age, by = "IID") %>% left_join(x = ., y = IPDGC_mds, by = c("FID","IID"))
  # 
  # 
  # NeuroX_mds <- fread(file = "NeuroX_MDS.mds", header = T) %>% select(., FID, IID, C1, C2, C3, C4)
  # NeuroX_age <- fread(file = "NeuroX_phenotype.txt", header = T, select = c(2,5)) %>% rename(., IID = SUBJECT_ID) 
  # 
  # IPDGC_cov <- left_join(x = FAM, y = IPDGC_age, by = "IID") %>% left_join(x = ., y = IPDGC_mds, by = c("FID","IID"))
  # NeuroX_cov <- left_join(x = FAM, y = NeuroX_age, by = "IID") %>% left_join(x = ., y = NeuroX_mds, by = c("FID","IID"))
  # 
  # fwrite(x = select(IPDGC_cov, -PID, -MID, -Sex, -Phenotype), file = "IPDGC.cov", row.names = F, sep = "\t")
  # fwrite(x = select(NeuroX_cov, -PID, -MID, -Sex, -Phenotype), file = "NeuroX.cov", row.names = F, sep = "\t")
}                                                      
library(tidyverse)
library(data.table)
# data pre-processing
{
  data_name <- c("IPDGC", "NeuroX", "PPMI")
  ExAC_maf <- fread(file = "/home/jinoo/skat-o/ExAC/ExAC_MAF.txt", header = T, data.table = F) %>% as_tibble() %>%
    rename(., ExAC_AF = AF_adj)
  test_fix <- list()
  
  # data select, ExAC, plink MAF
  {
    for(index in 1:3){
      {
      temp <- fread(file = paste0("/home/jinoo/skat-o/SKAT_data/",data_name[index], "_fix.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>%
        as_tibble() %>% arrange(CHROM) %>% mutate(., CPREA = paste0(CHROM,":",POS,REF,ALT))
      plink <- fread(file = paste0("/home/jinoo/skat-o/0311_test/",data_name[index],"_maf.frq"), header = T) %>%
        as_tibble() %>%select(., CHRPOS_plink = SNP, plink_maf = MAF)
      plink_case_control <- fread(file = paste0("/home/jinoo/skat-o/0311_test/",data_name[index],"_maf_case_control.frq.cc"), header = T) %>%
        as_tibble() %>%select(., CHRPOS_plink = SNP, plink_maf_case = MAF_A, plink_maf_control = MAF_U)


      CHRPOS2 <- list()
      test_id <- temp$ID;test_ch <- temp$CHROM;test_pos <- temp$POS
      for(change in 1:length(test_id)){
        if(test_id[change] == ""){
          CHRPOS2[[change]] <- paste0(test_ch[change], ":", test_pos[change])
        }
        else{
          CHRPOS2[[change]] <- test_id[change]
        }
      }
      test_fix[[index]] <- temp %>% bind_cols(., CHRPOS_plink = as.character(CHRPOS2)) %>%
        left_join(x = ., y = ExAC_maf, by = "CPREA") %>% left_join(x = ., y = plink, by = "CHRPOS_plink") %>%
        left_join(x = ., y = plink_case_control, by = "CHRPOS_plink")
      }
      
      # temp <- test_fix[[2]] %>% mutate(., CPREA = paste0(CHROM,":",POS,REF,ALT)) %>%
      #   left_join(x = ., y = ExAC_maf, by = "CPREA") %>%
      #   mutate(., ExAC_com = if_else(is.na(ExAC_AF), AF_adj, ExAC_AF)) %>%
      #   select(., -ExAC_AF, -AF_adj)

        
    }

  save(test_fix, file = "/home/jinoo/skat-o/SKAT_data/data_test_fix_0329.RData")
  }
}




data_name <- c("IPDGC", "NeuroX", "PPMI")
load(file = "/home/jinoo/skat-o/SKAT_data/data_test_fix_0329.RData")

#geneset load
{
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
  geneset_merge <- geneset_onehot %>% as.list()
  col_name <- names(geneset_merge)
}

#graph making!
for(index in 1:length(data_name)){
  # All ExAC compared each dataset
  {
    png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index], "_ExAC_compared.png"), width = 800, height = 600)
    plot(density(filter(test_fix[[index]], ExAC_AF >= 0.01 & ExAC_AF <= 0.03)$ExAC_AF, na.rm = T),
         main = paste0(data_name[index],"_MAF"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3, ylim = c(0,100))
    lines(density(filter(test_fix[[index]], plink_maf >= 0.01 & plink_maf <= 0.03)$plink_maf, na.rm = T), col = "blue", lty = 3, lwd = 3)
    legend("topright", c("ExAC", "plink"), fill = c("red","blue"), cex = 1.5)
    dev.off()
  }
  # ALL case-control
  {
    png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index], "_case_control_compared.png"), width = 800, height = 600)
    plot(density(filter(test_fix[[index]], plink_maf_case >= 0.01 & plink_maf_case <= 0.03)$plink_maf_case, na.rm = T),
         main = paste0(data_name[index],"_CASE_CONTROL_MAF(PLINK)"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3,ylim = c(0,100))
    lines(density(filter(test_fix[[index]], plink_maf_control >= 0.01 & plink_maf_control <= 0.03)$plink_maf_control, na.rm = T), col = "blue", lty = 3, lwd = 3)
    legend("topright", c("CASE", "CONTROL"), fill = c("red","blue"), cex = 1.5)
    dev.off()
  }
  
  # geneset_selection
  {
    for(gene in 2:length(col_name)){ #length(col_name)
      geneset <- geneset_merge[[gene]][!is.na(geneset_merge[[gene]])]
      print(paste0(data_name[index],"___",col_name[gene]))
      variant <- list()
      ## Parkinson geneset variant
      for(i in 1:length(geneset)){
        variant[[i]] <- subset(test_fix[[index]], subset = ((Gene.knownGene == geneset[i] & Func.knownGene == "exonic")))
      }
      variant <- bind_rows(variant)
      
      ### 1. nonsynonymous geneset
      {
        nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"
                                                   | ExonicFunc.knownGene == "stopgain"
                                                   | ExonicFunc.knownGene ==  "stoploss"
                                                   | ExonicFunc.knownGene ==  "frameshift_deletion"
                                                   | ExonicFunc.knownGene ==  "frameshift_insertion"
                                                   | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                                   | ExonicFunc.knownGene ==  "splicing"))
        
        # ExAC compared
        {
          png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index], "_",col_name[gene],"_ExAC_compared_nonsyn.png"), width = 800, height = 600)
          plot(density(filter(nonsynonymous, ExAC_AF >= 0.01 & ExAC_AF <= 0.03)$ExAC_AF, na.rm = T),
               main = paste0(col_name[gene],"_",data_name[index],"_Non-syn_MAF"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3, ylim = c(0,100))
          lines(density(filter(nonsynonymous, plink_maf >= 0.01 & plink_maf <= 0.03)$plink_maf, na.rm = T), col = "blue", lty = 3, lwd = 3)
          legend("topright", c("ExAC", "plink"), fill = c("red","blue"), cex = 1.5)
          dev.off()
        }
        
        # case-control
        {
          png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index], "_",col_name[gene], "_case-control_nonsyn.png"), width = 800, height = 600)
          plot(density(filter(nonsynonymous, plink_maf_case >= 0.01 & plink_maf_case <= 0.03)$plink_maf_case, na.rm = T),
               main = paste0(col_name[gene],"_",data_name[index],"_Non-syn_CASE_CONTROL_MAF(PLINK)"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3,ylim = c(0,100))
          lines(density(filter(nonsynonymous, plink_maf_control >= 0.01 & plink_maf_control <= 0.03)$plink_maf_control, na.rm = T), col = "blue", lty = 3, lwd = 3)
          legend("topright", c("CASE", "CONTROL"), fill = c("red","blue"), cex = 1.5)
          dev.off()
        }
        
      }  
      
      ### 2. CADD > 12.37 variant
      {
        cadd <- subset(variant, subset = ( ExonicFunc.knownGene == "nonsynonymous_SNV"
                                           | ExonicFunc.knownGene == "stopgain"
                                           | ExonicFunc.knownGene ==  "stoploss"
                                           | ExonicFunc.knownGene ==  "frameshift_deletion"
                                           | ExonicFunc.knownGene ==  "frameshift_insertion"
                                           | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                           | ExonicFunc.knownGene ==  "splicing") & CADD13_PHRED > 12.37)
        # ExAC compared
        {
          png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index], "_",col_name[gene], "_ExAC_compared_cadd.png"), width = 800, height = 600)
          plot(density(filter(cadd, ExAC_AF >= 0.01 & ExAC_AF <= 0.03)$ExAC_AF, na.rm = T),
               main = paste0(col_name[gene],"_",data_name[index],"_CADD_MAF"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3, ylim = c(0,100))
          lines(density(filter(cadd, plink_maf >= 0.01 & plink_maf <= 0.03)$plink_maf, na.rm = T), col = "blue", lty = 3, lwd = 3)
          legend("topright", c("ExAC", "plink"), fill = c("red","blue"), cex = 1.5)
          dev.off()
        }
      
        # case-control
        {
          png(filename = paste0("/home/jinoo/skat-o/0311_test/graph_maf/", data_name[index],"/",data_name[index],"_",col_name[gene],"_case-control_cadd.png"), width = 800, height = 600)
          plot(density(filter(cadd, plink_maf_case >= 0.01 & plink_maf_case <= 0.03)$plink_maf_case, na.rm = T),
               main = paste0(col_name[gene],"_",data_name[index],"_CADD_CASE_CONTROL_MAF(PLINK)"), cex.main = 1.5 ,col = "red", xlab = "MAF", lwd = 3, ylim = c(0,100))
          lines(density(filter(cadd, plink_maf_control >= 0.01 & plink_maf_control <= 0.03)$plink_maf_control, na.rm = T), col = "blue", lty = 3, lwd = 3)
          legend("topright", c("CASE", "CONTROL"), fill = c("red","blue"), cex = 1.5)
          dev.off()
        }
      }
      
        
        
    }
  }
  
  print(paste(data_name[index], "done!!!!", sep = " "))
}







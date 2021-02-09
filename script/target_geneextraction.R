geneset_extract_target <- function(target_gene, fix, type = "NonSyn"){
  print(paste0(target_gene,"_SetID making!!"))
  
  if(type == "NonSyn"){
    ## Parkinson geneset variant
    variant <- subset(fix, subset = ((Gene.knownGene == target_gene & Func.knownGene == "exonic")),
                               select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","PHRED","ID2"))
      
    ### 1. nonsynonymous geneset
    nonsynonymous <- subset(variant, subset = (ExonicFunc.knownGene == "nonsynonymous_SNV"))
    fwrite(x = subset(nonsynonymous, select = "ID2")[,1], file = paste0(target_gene,"_non-syn.txt"), col.names = F)
    
    ### 2. Lof (stopgain, stoploss, frameshift_deletion, frameshift_insertion, splicing, )
    lof <- subset(variant, subset = ( ExonicFunc.knownGene == "stopgain"
                                        | ExonicFunc.knownGene ==  "stoploss"
                                        | ExonicFunc.knownGene ==  "frameshift_deletion"
                                        | ExonicFunc.knownGene ==  "frameshift_insertion"
                                        | ExonicFunc.knownGene ==  "frameshift_block_substitution"
                                        | ExonicFunc.knownGene ==  "splicing"))
    fwrite(x = subset(lof, select = "ID2")[,1], file = paste0(target_gene,"_lof.txt"), col.names = F)
    
  } else { # ALL
    variant <- subset(fix, subset = ((Gene.knownGene == target_gene)),
                      select = c("CHROM", "POS", "ID", "REF","ALT", "Gene.knownGene","ExonicFunc.knownGene","PHRED","ID2"))
    fwrite(x = subset(variant, select = "ID2")[,1], file = paste0(target_gene,"_all.txt"), col.names = F)
  }
} 

plink_path <- "/home/dblab/temp/case_only_snptable/data/exclude_case"
target_gene <- "PRKN"
fix <- fix_load_KFPD(path = "/home/dblab/temp/case_only_snptable/data/exclude_case_fix_anno.txt",
                     plink_path = plink_path)
geneset_extract_target(target_gene = target_gene, fix = fix, type = "ALL")

system(glue("/tools/plink --bfile {plink_path} --extract {target_gene}_all.txt --allow-no-sex --recode A --out test_dosage_ALL",
            data_name = data_name, target_gene = target_gene, plink_path = plink_path), ignore.stdout = T)

all_dosage <- tibble()
tryCatch(expr = all_dosage <- fread(file = "test_dosage_ALL.raw", header = T) %>% select(-FID, -PAT, -MAT, -SEX) %>% as_tibble(),
         error = function(e){ print(e) })

View(all_dosage)
apply(X = all_dosage %>% select(-IID, -PHENOTYPE), MARGIN = 1, sum, na.rm = T) %>% tibble(count_sum = .) %>% 
  bind_cols(all_dosage %>% select(IID), .) %>% View()

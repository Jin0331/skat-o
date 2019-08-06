install.packages("DBI")
install.packages("RMySQL")

library(RMySQL)

driver <- dbDriver("MySQL")
con <- dbConnect(drv = driver, host = "210.115.229.96",dbname = "annotation",
                 user = "jinoo", pass = "sempre813!")

# query <- "select * from clinvar_20190609;"
# temp <- dbGetQuery(con, statement = query)

temp <- tbl(con, "cadd13")
chr17 <- temp %>% filter(CHROM == "17") %>%
  as_tibble

temp <- fread("whole_genome_SNVs.tsv.gz", sep = "\t", skip = 1)

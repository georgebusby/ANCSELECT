#######################################################
###  SCRIPT TO COMPARE PLS WITH GWAS SUMMARY STATS  ###
#######################################################

library( RSQLite )
library( hapdb )

# ## I RAN THIS FIRST BIT ON THE CLOSETS
summstatsfile <- "/mnt/kwiat/data/1/galton/malariagen/human/data/release/staging/public/EGAS00001001311/unencrypted/EGAS00001001311_MalariaGEN_GWAS_summary_statistics/EGAS00001001311_MalariaGEN_GWAS_summary_statistics.csv.gz"
#db <-  dbConnect( dbDriver( "SQLite" ), summstatsfile )
# summstats <- dbGetQuery( db, "SELECT `chromosome`, `position`, `rsid`, `Gambia_impute_info`, `Gambia_controls_alleleB_frequency`, `Kenya_impute_info`, `Kenya_controls_alleleB_frequency` FROM summary_statistics WHERE chromosome != '0X'" )
# write.table(summstats, file = 'data/KenyaGambiaAlleleFreqs.txt.gz', col.names = T, row.names = F)
summstats <- read.csv(summstatsfile,header = T)


## GET IHH DATA AND COMPARE TO UIHS
ihh <- read.table('~repos/glycophorings/data/KenyaDUP4QuangIHS.txt.gz', header = T)
pls <- read.table()





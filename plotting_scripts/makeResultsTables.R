#######################################################
### SCRIPT TO PLOT LRT ACROSS GWAS ETHNIC GROUPS    ###
#######################################################
library("dplyr")
source("plotting_scripts/gavinManhattenFunctions.R")

#############################################
## GET DATA
snps <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330Kphased.legend", header=F,as.is=T)
colnames(snps) <- c("chrom","rsid","pos","a0","a1")
popkey_file <- "data/PopulationKey.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$AncestryRegion <- gsub("Afro-Asiatic","Afroasiatic",popkey$AncestryRegion)
colour_table <- read.table("data/RegionColours.txt", header = T, stringsAsFactors = F, comment.char = "")
ancreg_list <- colour_table$AncestryRegion
pcolshex <- colour_table$Colour
#############################################
## GET LRT RESULTS
pop <- "FULAI"
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                            "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                            "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

poplist <- popkey$Ethnic_Group[1:48]


allhits <- c()
for(pop in poplist)
{
 
  hit_tab <- paste0("data/",pop,"LRTmanhattanIII_hits.txt")
  ginf_tab <- paste0("data/",pop,"LRTmanhattanIII_ginfo.txt")

  if(file.exists(hit_tab))
  {
    hits <- read.table(hit_tab, sep = "\t", header = T, stringsAsFactors = F, as.is = T)
    hits <- cbind(pop,hits,0)
    colnames(hits)[1] <- "pop"
    colnames(hits)[ncol(hits)] <- "genes"
    ginfo <- read.table(ginf_tab, sep = "\t", header = T, stringsAsFactors = F, as.is = T)  
    
    ## 
    for(i in hits$rsid)
    {
      tmpgenes <- ginfo$name2[which(ginfo$rsid == i)]
      if(is.na(tmpgenes)) tmpgenes <- "No nearby genes"
      if(length(tmpgenes) > 1) tmpgenes <- paste(tmpgenes, collapse = ";")
      hits$genes[which(hits$rsid==i)] <- tmpgenes
    }
    allhits <- rbind(allhits,hits)
  }
}


## NUMBER OF POPS WITH HITS > 8
sum(sort(table(allhits$pop[allhits$P>=8 & allhits$keep == 1])) > 0)
sum(sort(table(allhits$pop[allhits$P>=8 & allhits$keep == 1])) > 5)

test <- table(allhits$pop[allhits$P>=8 & allhits$keep == 1]) >0 & table(allhits$pop[allhits$P>=8 & allhits$keep == 1]) < 5

pops2 <- names(table(allhits$pop[allhits$P>=8 & allhits$keep == 1])[test])

allhits[which(allhits$pop %in% pops2),1:8]


### MAKE NICE PRETTY TABLE FOR SUPPLEMENT ##

allhits_print <- allhits[allhits$keep == 1,c(1,2,3,4,5,6,9,10,11,13)]
allhits_print$generegion <- paste0(allhits_print$chrom,": ",
                                   round(allhits_print$lowerpos/1e6,2),"-",
                                   round(allhits_print$upperpos/1e6,2),"Mb")
allhits_print$region <- gsub("_"," ",allhits_print$region)
allhits_print <- allhits_print[,c(1,11,3,5,6,7,10)]
allhits_print$pop <- tidyNames(as.character(allhits_print$pop), fula = T)
allhits_print$genes <- gsub(";","; ",allhits_print$genes)
allhits_print <- xtable(allhits_print,align="|r|l|c|c|c|c|l|l|", digits = 3)
### print all events
newlines <- c()
newlineord <- allhits_print$pop
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "main results table"
addtorow$command[[2]] <- "\\hline \n"
print(allhits_print, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=FALSE,
      include.rownames=FALSE,
      caption.placement="top",file="data/AllPopsSignifHits.tex",
      booktabs=TRUE,add.to.row=addtorow)#,



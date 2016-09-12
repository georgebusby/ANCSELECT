#####
library(dplyr) ## for binding dataframes with difference numbers of columns; bind_rows

in_dir <- '/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/'
chroms <- c(1:22)
pops <- read.table('/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP.idfile.txt')
pops <- as.character(levels(pops[,2]))
regions <- read.table('/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/FULAInolocal.idfile2.txt') 
regions <- as.character(levels(regions[,2]))
regions <- sort(c("East_Africa_Afroasiatic",regions))

tophits <- c()
pops<- "FULAI"
for(pop in pops)
{
  print(pop)
  liks <- c()
  for(chrom in chroms)
  {
    in_file <- paste0(in_dir,pop,'nolocalChrom',chrom,'PP.likelihoods.gz')
    tmp <- read.table(in_file,header=T,as.is=T)
    liks <- rbind(liks,tmp) 
  }
  liks <- cbind(pop,liks)

  
  ## GET TOP HITS
  pcols <- grep(".P",colnames(liks))
  tmp <- liks[rowSums(liks[,pcols]>6)>0,]
  tophits <- bind_rows(tophits,tmp)
}

## REMOVE EURASIAN DATA FOR THE TIME BEING
eurasia <- c('CDX','CEU','CHB','CHS','GBR', 'GIH','FIN','IBS','JPT','KHV','PELII','TSI')
tophits2 <- tophits[!tophits$pop%in%eurasia,]
## EXPLORE
table(tophits2$pop)
## only 20 groups have hits > 6

### FOR EACH POP, WORK OUT
## A: NUMBER OF DISTINCT SIGNALS
## B: ANCESTRY OF DISTINCT SIGNAL
## C: STRENGTH OF DISTINCT SIGNAL
## D: GENES IN DISCINT SIGNAL


source("~/R/Copy/Rprojects/AfricaPOPGEN/functions/hitplots.R")
genes <- load.genes()

genepops <- unique(tophits3$pop)
## REMOVE DUPLICATE SNPS FROM SAME GENETIC REGION
genelist <- c()
for(i in genepops)
{
  tmp <- tophits2[tophits2$pop==i,]
  for(chrom in unique(tmp$chrom))
  {
    tmp2 <- tmp[tmp$chrom==chrom ,]
    tmp2$pos2 <- round(tmp2$pos,-6)
    for(j in unique(tmp2$pos2))
    {
      tmp3 <- tmp2[tmp2$pos2==j,]
      ## FIND BIGGEST P VALUE IN REGION
      whichp <- apply(tmp3[,grep(".P",colnames(tmp3))],2, max)
      whichreg <- names(sort(whichp,decreasing=T)[1])
      whichp <- as.numeric(sort(whichp,decreasing=T)[1])
      tmp3 <- cbind(tmp3[tmp3[,whichreg] == whichp,][1,],gsub(".P","",whichreg))
      colnames(tmp3)[ncol(tmp3)] <- "ancestry.region"
      tmp3 <- cbind(tmp3,whichp)
      colnames(tmp3)[ncol(tmp3)] <- "ancestry.region.P"
      genelist <- rbind(genelist,tmp3)
    }
  }
}
   
genelist <- genelist[order(genelist$chrom,genelist$pos2),]
#genelist <- genelist[genelist$ancestry.region.P>=8,]

table(genelist$pop)

test <- genelist[genelist$pop=="WOLLOF",]
test <- test[order(test$ancestry.region.P,decreasing = T),]
head(test)


#####
## GENERALLY SPEAKING, APART FROM A COUPLE OF EXCEPTIONS, IT'S APPROPRTIAT
## TO CONCNTRATE ON FULA, JOLA AND JUHOAN AS THESE ARE THE THREE GROUPS WIHT
## MULTIPLE SIGNALS -LOG10 P VALUE > 8
## OTHER THINGS TO MENTION ARE:
##      LDB2 IN EAST AFRICA
##      HLA IN THE GAMBIA
##      NRG IN ARI?


## COMPARE TO RYAN'S CODE









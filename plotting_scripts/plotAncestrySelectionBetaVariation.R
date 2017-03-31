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

pop <- "FULAI"
chroms <- 1:22
# in_root <- paste0("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/ancestry_selection/",pop,"nolocalChrom")
in_root <- paste0("/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/",pop,"nolocalChrom")
in_suff <- ".ancestryselectionIII.gz"
pl <- c()
test <- file.exists(paste0(in_root,"01",in_suff))
if(test)
{
  for(chrom in chroms)
  {
    if(chrom < 10) chrom <- paste0("0",chrom)
    print(paste("loading chromosome:",chrom))
    in_file <- paste0(in_root,chrom,in_suff)
    if(chrom == 1)
    {
      pl <- read.csv(in_file,header=T, as.is=T)
    } else
    {
      pl <- rbind(pl,read.csv(in_file,header=T, as.is=T))
    }
  }
}

## COVARIANCE/CORRELATION OF BETAS ACROSS ANCESTRIES
self_reg <- popkey$AncestryRegion[popkey$Ethnic_Group==pop]
covmat <- matrix(0,nc = length(ancreg_list)-1,
                 nr = length(ancreg_list)-1)
rownames(covmat) <- colnames(covmat) <- ancreg_list[ancreg_list!=self_reg]
threshold <- 5
propthreshold <- 5
for(i in 1:nrow(covmat))
{
  ii <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.beta")]
  propi <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".prop")] >= propthreshold
  
  covmat[i,i] <- var(ii[propi], na.rm = T)
  for(j in 1:ncol(covmat))
  {
    if(i != j & j > i)
    {
      ## look at subset of SNPs where either ancestry is significant
      jj <- pl[,paste0(gsub("\\-","\\.",colnames(covmat)[j]),".RC.beta")]
    
      signi <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.P")] >= threshold
      signj <- pl[,paste0(gsub("\\-","\\.",colnames(covmat)[j]),".RC.P")] >= threshold
      propj <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".prop")] >= propthreshold
      
      
      s <- signi | signj
      s[is.na(s)] <- FALSE
      p <- propi & propj
      p[is.na(p)] <- FALSE
      
      ii[is.na(ii)] <- 0
      jj[is.na(jj)] <- 0
 
#      covmat[i,i] <- var(ii)
      covmat[i,j] <- cov(ii[which(s & p)],jj[which(s & p)])
      covmat[j,i] <- cov(ii[which(!s & p)],jj[which(!s & p)])
    }
  }
}
covmat <- covmat[,rev(colnames(covmat))]

##########################################################

png("figures/FULAILRTcovarianceDensities.png",
    height = 500, width = 1000)
layout(matrix(c(1,2,
                1,3),2,2, byrow = T))
## plot heat matp
par(mar = c(10,10,1,1))
image(1:nrow(covmat),
      1:ncol(covmat),
      covmat, axes = F,
      xlab = "", ylab = "")

for(i in 1:nrow(covmat))
{
  ii <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.beta")]
  for(j in 1:nrow(covmat))
  {
    
    if(rownames(covmat)[i] == colnames(covmat)[j])
    {
      text.lab <- paste0(signif(covmat[i,j],3),"\n[",signif(mean(ii,na.rm = T),3),"]")
    } else
    {
      text.lab <- signif(covmat[i,j],3)
    }
    text(i,j, labels = text.lab)
  }
}

axis(1,at = 1:nrow(covmat), labels = gsub("Africa_","Africa\n",rownames(covmat)), las = 2,lwd = 0)
axis(2,at = 1:nrow(covmat), labels = gsub("Africa_","Africa\n",colnames(covmat)), las = 2,lwd = 0)


  
## plot denties

par(mar = c(4,4,3,1))
plot(0,0,xlim= c(-2,2),ylim= c(0,1.6),xaxs = "i",yaxs = "i",
     type = "n",axes = F, xlab = "", ylab = "")
axis(1)
#axis(2, las = 2)
abline(v=0)

for(i in 1:nrow(covmat))
{
  ii <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.beta")]
  plotcol <- pcolshex[ancreg_list == rownames(covmat)[i]]
  s <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.P")] >= threshold
  p[is.na(p)] <- FALSE
  points(density(ii[p], na.rm = T)$x,
         density(ii[p], na.rm = T)$y,
         lwd = 2, col = plotcol, type = "l")
}

title(xlab = expression(italic(beta)))
title(main = "Genome-wide distribution")

## plot denties of significant betas

plot(0,0,xlim= c(-2,2),ylim= c(0,1.6),xaxs = "i",yaxs = "i",
     type = "n",axes = F, xlab = "", ylab = "")
axis(1)
#axis(2, las = 2)
abline(v=0)

for(i in 1:nrow(covmat))
{
  ii <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".RC.beta")]
  plotcol <- pcolshex[ancreg_list == rownames(covmat)[i]]
  p <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".prop")] >= propthreshold
  p <- pl[,paste0(gsub("\\-","\\.",rownames(covmat)[i]),".prop")] >= propthreshold
  p <- s & p
  p[is.na(p)] <- FALSE
  points(density(ii[p], na.rm = T)$x,
         density(ii[p], na.rm = T)$y,
         lwd = 2, col = plotcol, type = "l")
}

title(xlab = expression(italic(beta)))
title(main = "Significant SNPs")


dev.off()


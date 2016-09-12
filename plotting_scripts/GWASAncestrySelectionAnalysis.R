#######################################################
### SCRIPT TO WORK WITH COPYING PROBS FOR GWAS INDS ###
#######################################################

## THE FILES ARE STORED AS HDF5 FILES ON /well AND THIS
## SCRIPT HELPS TO PULL OUT THE RELEVANT INFO

library(h5) ## INSTALL IF NECESSARY
#library(abind) ## for combining arrays
## WE'LL WORK ON A SINGLE CHROMOSOME, COPYPROBS ARE SPLIT ACROSS
## 22 FILES, ONE PER CHROMOSOME
chrom <- "02"
in_dir <- "/mnt/kwiat/well/human/george/GWAS_painting/chromopainter/hdf5/"
## in_dir IS THE RELEVANT DIRECTORY MOUNTED ON MY OWN COMPUTER -- THIS FULL PATH IS BELOW
#in_dir <- "/well/malariagen/malariagen/human/george/GWAS_painting/chromopainter/hdf5/"
hdf5file <- paste0(in_dir,"GambiaChrom",chrom,"Copyprobsperlocus.h5")

## DEFINE OUR DATASET
datafile <- h5file(hdf5file, mode = "r")

## THE BASIC LAYOUT OF THIS FILE IS:
##    /copyprobsperlocus/chrom04/
##    /copyprobsperlocus/chrom04/donors  = painting donors
##    /copyprobsperlocus/chrom04/map  = info on snps
## THEN WE HAVE THE ACTUAL DATA STORED IN CHUNKS OF 200 HAPS
## TO GET AT THE INDS, WE LOOK HERE
##    /copyprobsperlocus/chrom04/haps/1-100 etc
## THE ACTUAL COPYIN PROBS ARE HERE
##    /copyprobsperlocus/chrom04/probs/1-100 etc


## LET'S TRY WITH A SINGLE CHROMOSOME: 2
## WE NEED:
## AVERAGE GENOME-WIDE COPYING PER HAPLOTYPE

#######################################################
## 00 THE FIRST THING IS TO GET THE SNP INFO
## GET MAP AND POSITION INFO
snps <- data.frame(readDataSet(datafile[paste0("/copyprobsperlocus/chrom",chrom,"/map")]), stringsAsFactors = FALSE)
colnames(snps) <- c("rsid","recrateCP","chrom","position", "a0","a1")
snps$position <- as.numeric(snps$position)
## SELF EXPLANATORY, I HOPE -- recrateCP IS THE RECOMBINATION RATE BETWEEN
## A SNP AND THE NEXT ONE, MEASURED IN MORGANS AND USED BY CHROMOAPAINTER

#######################################################
## 01 LET'S WORK OUT WHICH INDS ARE IN THIS FILE
inds_available <- list.datasets(datafile[paste0("/copyprobsperlocus/chrom",chrom,"/haps")])

## MAKE A DATAFRAME OF OUR DATASETS: THIS RELATES INDIVIDUAL IDS TO THEIR PLACE IN OUR DATASET
haps <- c()
for(i in inds_available)
{
  tmp <- data.frame(readDataSet(datafile[i]))
  inds <- strsplit(i, split="\\/")[[1]][5]
  ind1 <- as.numeric(strsplit(inds, split="-")[[1]][1])
  ind2 <- ind1 + (nrow(tmp)-1)
  haps <- rbind(haps,cbind(i,ind1:ind2,1:nrow(tmp),tmp))
}
colnames(haps) <- c("hdf5dataset","hdf5_ID", "hdf5dataset_ID","haplotypeID")
head(haps)
## WHEN I RUN THIS, I GET
# hdf5dataset hdf5_ID hdf5dataset_ID             haplotypeID
# 1 /copyprobsperlocus/chrom04/haps/1-100       1              1 HAP 1 6006278010_R01C01
# 2 /copyprobsperlocus/chrom04/haps/1-100       2              2 HAP 2 6006278010_R01C01
# 3 /copyprobsperlocus/chrom04/haps/1-100       3              3 HAP 1 6006278010_R02C01
# 4 /copyprobsperlocus/chrom04/haps/1-100       4              4 HAP 2 6006278010_R02C01
# 5 /copyprobsperlocus/chrom04/haps/1-100       5              5 HAP 1 6006278010_R03C01
# 6 /copyprobsperlocus/chrom04/haps/1-100       6              6 HAP 2 6006278010_R03C01
## WHICH TELLS ME THAT HAPLOTYPE 1 IS HAP1 OF IND 6006278010_R01C01 AND IS STORED A POSITION 1 IN hdfdataset
##                               2 IS HAP2 OF IND 6006278010_R01C01 AT POSITION 2, etc
## LET'S GET THE CASE/CONTROL STATUS
casecontrol <- read.table("data/GWAScasecontrol.txt",header = T, as.is = T)
casecontrol <- casecontrol[-1,]
hap_casecontrol <- casecontrol[match(sapply(as.character(haps$haplotypeID), function(x){strsplit(x, split=" ")[[1]][3]}),casecontrol[,1]),2]

## LET'S GET ETHNICITY
ethnicity  <- read.table("data/GWASethnicity.txt", header = T, as.is = T)
ethnicity <- ethnicity[-1,]
hap_ethnicity <- ethnicity[match(sapply(as.character(haps$haplotypeID), function(x){strsplit(x, split=" ")[[1]][3]}),ethnicity[,1]),4]
hap_ethnicity[is.na(hap_ethnicity)] <- "OTHER"

ethnic <- "FULA"
ethnic_indices <- which(hap_ethnicity==ethnic)
#######################################################
## 01 LOOK AT A SNPS ACROSS DIFFERENT INDS FROM ETHNIC 
snp_index <- 1:nrow(snps)
## TO PULL OUT THE COPYING PROBS FOR THIS REGION, WE DO THE FOLLOWING
## THE PROBS ARE STORED AS A 3D ARRAY: n_haps * n_snps * n_dons
probs_available <- list.datasets(datafile[paste0("/copyprobsperlocus/chrom",chrom,"/probs/")])

## WE CAN GET THE DONORS HERE
donors <- datafile[paste0("/copyprobsperlocus/chrom",chrom,"/donors")][]
## WE CAN CONVERT THE DONORS TO ANCESTRY REGIONS LIKE THIS
popkey <- read.table("data/PopulationKey.txt", header = T, as.is = T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
donor_regions <- popkey$AncestryRegion[match(donors,popkey$Ethnic_Group)]
regions <- sort(unique(popkey$AncestryRegion))

#######################################################
## 02 WHAT ARE THE AVERAGE COPYING PROBS FROM EACH ANCESTRY REGION
##    FOR EACH ETHNIC HAPLOTYPE

ethnic_haps <- haps[ethnic_indices,]
probs_haps <- matrix(0,nr=0,nc=length(donors))
probs_snps <- matrix(0,nr=nrow(snps),nc=length(donors))
for(i in as.character(unique(ethnic_haps$hdf5dataset)))
{
  hap_index <- ethnic_haps$hdf5dataset_ID[ethnic_haps$hdf5dataset==i]
  tmp <- datafile[gsub("haps","probs",i)][hap_index,,]
  
  ## GET CHROMOSOME-WIDE COPYING PROBS
  probs_tmp <- apply(tmp,c(1,3),mean)
  probs_tmp <- probs_tmp/rowSums(probs_tmp)
  probs_haps <- rbind(probs_haps,probs_tmp)

  ## GET COPYING PROBS AT SNPS
  probs_tmp <- apply(tmp,c(2,3),sum)
  probs_tmp <- (probs_tmp/rowSums(probs_tmp))*dim(tmp)[1]
  probs_snps <- probs_snps + probs_tmp
  
  print(i)  
}

###################################################################
## WRITE TO FILE TO SAVE ##
# write.table(probs_snps, file = 'data/FULACopyingProbsPerSNPChromosome02.txt',
#             quote = F, row.names = F, col.names = F)
# write.table(probs_haps, file = 'data/FULACopyingProbsChromosome02.txt',
#             quote = F, row.names = F, col.names = F)
###################################################################
# prob_haps <- read.table('data/FULACopyingProbsPerSNPChromosome02.txt', header = F, as.is = T)
# prob_snps <- read.table('data/FULACopyingProbsChromosome02.txt', header = F, as.is = T)
###################################################################

## NOW SUM ACROSS COLUMNS (IE DONORS FROM THE SAME REGION)
probs_haps_reg <- array(0, dim = c(nrow(probs_haps),length(regions)))
for(i in regions)
{
  if(i %in% unique(donor_regions))
  {
    probs_haps_reg[,regions==i] <- apply(probs_haps[,donor_regions==i],1,sum)
  }
}
colnames(probs_haps_reg) <- regions

## SUM ACROSS COLUMNS AT SNPS
probs_snps_reg <- array(0,dim = c(nrow(probs_snps),length(regions)))
for(i in regions)
{
  if(i %in% unique(donor_regions))
  {
    probs_snps_reg[,regions==i] <- apply(probs_snps[,donor_regions==i],1,sum)
  }
}
colnames(probs_snps_reg) <- regions
probs_snps_reg <- probs_snps_reg/rowSums(probs_snps_reg)

#################################################
## LET'S LOOK AT A SECTION OF CHROMOSOME 2 TO SEE
## IF THERE ARE LARGE DIFFERENCES BETWEEN INDS GENOME-WIDE
## COPYING AND THAT AT A SNP
tvd <- function(x,y){
  ## x is a vector of probs, e.g. a hap's copying probs at a snp
  ## y is another vector of probs, e.g. a hap's genome-wide copying probs
  tvd <- 0.5*sum(abs(x - y))
  return(tvd)
}

## SNP INDEX
snp_index <- which(snps$position>=130e6&snps$position<=140e6)

## A MATRIX OF L SNPS * N HAPS
probs_tvd <- matrix(0,nr=length(snp_index),nc=nrow(probs_haps))
cnt <- 1
for(i in as.character(unique(ethnic_haps$hdf5dataset)))
{
  hap_index <- ethnic_haps$hdf5dataset_ID[ethnic_haps$hdf5dataset==i]
  tmp <- datafile[gsub("haps","probs",i)][hap_index,,]
  tmp <- tmp[,snp_index,]
  
  hap_probs <- probs_haps_reg[hap_index,]
  
  for(j in 1:nrow(hap_probs))
  {
    tmp2 <- tmp[j,,]
    tmp3 <- matrix(0,nr=nrow(tmp2),nc=length(regions))
    ## CONVERT TO DONOR REGIONS
    for(k in regions)
    {
      if(k %in% unique(regions))
      {
        tmp3[,regions==k] <- apply(tmp2[,donor_regions==k],1,sum)
      }
    }
    colnames(tmp3) <- unique(regions)
    tmp3 <- tmp3[,colnames(hap_probs)]
    tvds <- 0.5*apply(apply(tmp3,1,function(x){abs(x-hap_probs[j,])}),2,sum)
    probs_tvd[,cnt] <- tvds
    cnt <- cnt + 1
  }
}
  
meds <- apply(probs_tvd,1,median)
int1 <- apply(probs_tvd,1,quantile,0.25)
int2 <- apply(probs_tvd,1,quantile,0.75)

xs <- snps$position[snp_index]
plot(xs,meds,ylim=c(min(int1),max(int2)), axes = F, xlab = "position on chromosome 2",
     ylab = "TVD", type = "n")
axis(2, las = 2)
axis(1)
polygon(c(xs,rev(xs)), c(int1,rev(int2)), col = "grey", border = NA)
points(xs,meds, col = "black", type = "S", lwd = 1)
#points(xs,apply(probs_tvd,1,max), col = "red", type = "S", lwd = 2)
#points(xs,apply(probs_tvd,1,min), col = "red", type = "S", lwd = 2)
box()



  
  ## GET CHROMOSOME-WIDE COPYING PROBS
  probs_tmp <- apply(tmp,c(1,3),mean)
  probs_tmp <- probs_tmp/rowSums(probs_tmp)
  probs_haps <- rbind(probs_haps,probs_tmp)
  
  ## GET COPYING PROBS AT SNPS
  probs_tmp <- apply(tmp,c(2,3),sum)
  probs_tmp <- (probs_tmp/rowSums(probs_tmp))*dim(tmp)[1]
  probs_snps <- probs_snps + probs_tmp
  
  print(i)  
}


#################################################
## NOW PLOT TO SEE WHAT IT LOOKS LIKE
## GET MY USUAL COLOUR SCHEME
colour_table <- read.table("data/RegionsColours.txt", header = F, as.is = T, comment.char="")
region_order <- colour_table[,2]
region_colours <- colour_table[,1]

png("figures/Chromosome02copyingProbs.png",width = 2000, height = 500, res = 150)
  par(mar=c(4,2,1,1))
  barplot(t(probs_snps_reg),border=NA,space=0,
          width=diff(snps$position),
          col=region_colours[match(region_order,colnames(probs_snps_reg))],
          axes = F, xlab="position on chromosome (Mb)")
  
  x_at <- seq(0,max(snps$position),25e6)
  if(max(x_at) < max(snps$position))
  {
    x_at <- c(x_at,max(snps$position))
  }else
  {
    x_at[length(x_at)] <- max(snps$position)
  }
  axis(1, at = x_at, labels = round(x_at/1e6))
dev.off()


####################
## RYAN'S CODE ...
## LOGISTIC FUNCTIONS
####### NEW STANDARD, COOKBOOK APPROACH
## 00 DEFINE A MATRIX OF MUS
# mus <- probs_haps_reg
# ## remove column of 0s
# mus <- mus[,apply(mus,2,mean)!=0]
# ## LOGIT OF MUS
# mus <- log(mus/( 1-mus))
# find_beta <- function(beta,sum_y,mu){return(sum_y - sum( 1/(1+exp(-(mu+beta)))))}
# sigma_sq<-2
# 
# out_mat <- snps
# for(donor_reg in colnames(mus))
# {
#   print(donor_reg)
#   mu_i <- mus[,donor_reg] 
#   ps <- lrs <- rep(0,nrow(out_mat))
#   for(snp_index in 1:nrow(snps))
#   {
#     sum_yi <- probs_snps_reg[snp_index,donor_reg]
#     beta_hat<-uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
#     LRT<-2*(beta_hat*sum_yi+sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat)))))
#     lrs[snp_index] <- LRT
#     ps[snp_index]<- -log10(pchisq(LRT,df = 1,lower.tail = F, log.p = T))
#   }
#   out_mat <- cbind(out_mat,lrs,ps)
#   colnames(out_mat)[ncol(out_mat)-1] <- paste(donor_reg,"likelihood",sep=".")
#   colnames(out_mat)[ncol(out_mat)] <- paste(donor_reg,"P",sep=".")
# }
  
  







## NOW GET 



################################################################
## SPLIT BY CASE/CONTROLS ANF POPULATION AND THEN ORDER HAPS
## ORDER HAPS BY MOST VARIABLE ANCESTRY REGION TO MAKE PLOT PRETTY
pdf("figures/GambiaGWASChromosome2LCTCopyingProbs.pdf", height = 7, width = 10)
layout(matrix(c(1:13,13),7,2, byrow=T))
for(cc in c("CONTROL", "CASE"))
{
  plot_mat <- new_probs[hap_casecontrol == cc,region_order]
  hap_order <- order(plot_mat[,which.max(apply(plot_mat,2,var))])
  plot_mat <- plot_mat[hap_order,region_order]
  plot_cols <- region_colours
  par(mar=c(0.5,1,3,1))
  bp <- barplot(t(as.matrix(plot_mat)),
                yaxt="n",xlab="",beside=F,
                main=paste(cc, "[",nrow(plot_mat)/2, "inds]"),cex.main=0.75,border=NA,horiz=F,
                col=plot_cols,xaxt="n",add=F,plot=T)
}
big_ethnicities <- names(which(table(hap_ethnicity)>100))
for(cc in big_ethnicities)  
{
  plot_mat <- new_probs[hap_ethnicity == cc,region_order]
  eth_casecontrol <- casecontrol[match(sapply(as.character(haps$haplotypeID[hap_ethnicity==cc]), function(x){strsplit(x, split=" ")[[1]][3]}),casecontrol[,1]),2]
  hap_order <- order(eth_casecontrol,plot_mat[,which.max(apply(plot_mat,2,var))], decreasing = T)
  plot_mat <- plot_mat[hap_order,region_order]
  plot_cols <- region_colours
  par(mar=c(0.5,1,3,1))
  bp <- barplot(t(as.matrix(plot_mat)),
                yaxt="n",xlab="",beside=F,
                main=paste(cc, "[",nrow(plot_mat)/2, "inds]"),cex.main=0.75,border=NA,horiz=F,
                col=plot_cols,xaxt="n",add=F,plot=T)
  ## WORK OUT CASECONTROL
  for(i in c("CONTROL","CASE"))
  {
    axcol <- "black"
    if(i == "CASE") axcol <- "red"
    xat <- which(eth_casecontrol[hap_order]==i)
    axis(1,at=bp[xat],labels=NA, col.ticks=axcol)
    
  }    
}

## LEGEND
plot(0,0,type="n", xlab="n",ylab="", axes = F)
legend("left", legend=gsub("\\_"," ",region_order), fill=region_colours, border = region_colours, ncol = 3, bty="n")
legend("bottomright",legend = c("CONTROL", "CASE"), lty = 1, col = c("black","red"), bty="n")

dev.off()
## NO REAL CLEAR TRENDS FROM THE PLOTS THAT ONE ONE COPYING IS ASSOCIATED WITH CASES OR CONTROLS ##



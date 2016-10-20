###
### SCRIPT TO PLOT SHARED SIGNAL OF ANCESTRY DEVIATION AT HLA IN THE GAMBIA ###
#setwd(paste0(main_dir,"ANCSELECT/"))
library("h5")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"ANCSELECT/"))
###########################################################
## DEFINE DATAFILES
#fsanalyname <- 'MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalForce.A'
#leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
#poppos_file <- "data/MalariaGenAdmixturePopulationKeyMapPositions.txt"
#popkey_file <- "data/PopulationKey.txt"
#popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
#pca_file <- "data/Africa300KPCS.txt"
#tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
#mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")
datafile <- h5file('/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5')


## LOAD POPKEY FILE ##
popkey_file <- "~/repos/admixture_in_africa/data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                            "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                            "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

analysis <- "nonlocal"


#####
## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
  ## finds the first sample index for a haplotype
  sample <- (hap*nsamps)-(nsamps-1)
  return(sample)
}


pop <- "FULAI"
chrom <- "02"
plot_range <- c(130e6,140e6)
gene_region <- c(136e6,137e6)


pop="JOLA"
chrom <- "06"
plot_range <- c(30e6,35e6)
gene_region <- c(32.2e6,33e6)

## GET MAP AND POSITION INFO
map <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/map")]),stringsAsFactors = F)
colnames(map) <- c("position","recrate")
map <- data.frame(apply(map,2,as.numeric))
snps <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/snps")]))
colnames(snps) <- c("chrom","rsid","pos","a0","a1")
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- readDataSet(datafile[paste0("/paintings/samples/individuals")])
colnames(psamples) <- c("ind","region","X")

####################################################################################
## GET PAINTINGS ACROSS ALL POP OF INTEREST
psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
# 2 haps per sample!!
psampleshap <- hap2sampleindex(psamplesind,2)
psampleshap <- sort(c(psampleshap,psampleshap+1))
psamplesindsamp <- hap2sampleindex(psampleshap)
tmp <- c()
for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
paintedchrom <- datafile[paste0("/paintings/chrom",chrom,"/",analysis)]
paintedchrom <- paintedchrom[,tmp]

## SWITCH DONORS TO REGIONS
happops <- c()
for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
hapregs <- c()
for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
hapregs <- factor(hapregs,levels=ancreg_list)
paintedchromreg <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchromreg[,i] <- hapregs[paintedchrom[,i]]
paintedchrompop <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchrompop[,i] <- happops[paintedchrom[,i]]

## GET COPIER COPIES ACROSS CHROMOSOME
copiercopies <- datafile[paste0("/lengths/chrom",chrom,"/copiercopies")]
# copiercopiesreg <- matrix(nc=ncol(copiercopies),nr=nrow(copiercopies))
# for(i in 1:ncol(copiercopies)) copiercopiesreg[,i] <- hapregs[copiercopies[,i]]
# copiercopiespop <- matrix(nc=ncol(copiercopies),nr=nrow(paintedchrom))
# for(i in 1:ncol(copiercopies)) copiercopiespop[,i] <- happops[copiercopies[,i]]

# ## FOR EACH VALUE IN PAINTEDCHROM, WE WANT TO 
# ## GET THE ID OF THE PERSON WHO COPIER COPIES
# pchrom <- paintedchrom[plot_region,]
# for(i in 1:nrow(pchrom))
# {
#   j <- which(plot_region)[i]
#   print(i)
#   pchrom[i,] <- copiercopies[j,pchrom[i,]]
# }

## NOW RUN TEST AT EACH SNP MASKING THE DATA

loglik <- function(v,data,nsamps=n_samps,useall=useallsamps,slike=samp2use)
{
  ## GENERATE AN IDENTITY MATRIX WHERE 1ST COL IS 1:nhaps
  ## AND NEXT nsamps COLS ARE THE IDENTITY OF THE DONOR THAT
  ## EACH HAPLOTYPE COPIES AT THIS SNP
  nhaps <- length(data)/nsamps
  ind <- matrix(cbind(1:nhaps,matrix(data,ncol=nsamps,byrow=T)),ncol=nsamps+1)
  ## MAKE AN EMPTY MATRIX AND FILL WITH THE PROBABILITIES
  ## OF COPYING FROM THAT REGION AT THAT SNP ie BY USING v
  prob <- matrix(NA, nrow=nhaps,ncol=nsamps)
  for(k in 2:(nsamps+1)) prob[,(k-1)] <- v[ind[,c(1,k)]]
  
  ## NOW GENERATE LOG LIKELIHOOD ACROSS ALL HAPLOTYPES
  ## CURRENT IMPLEMENTATION USES A SINGLE PAINTING: slike DEFINES THE SAMPLE TO USE
  if(useall == FALSE ) prob <- sum(log(prob[,slike]))
  #if(useall == TRUE ) prob <- sum(log(rowSums(prob)/ncol(prob)))
  ## THIS IS GBJB'S PREVIOUS VERSION, WHICH LOOKS LIKE IT WORKS
  ## BUT IS PERHAPS INCORRECT THING TO DO
  if(useall == TRUE) prob <- mean(apply(prob,2,function(x)sum(log(x), na.rm = T)), na.rm = T)
  return(prob)
}

par.loglik <- function(v,data,nsamps=n_samps,useall=useallsamps,slike=samp2use,lambda,colindex)
{
  ## WE NEED TO CHOOSE A REGION TO PIVOT THE LAMBDA VALUE AROUND
  ## WE WILL USE COLUMN 1, UNLESS THIS IS THE REGION THAT THE POP
  ## OF INTEREST COMES FROM, IN WHICH CASE THE COPYING PROPS WILL
  ## BE 0 (AS WE DIS-ALLOW) SELF-COPYING. IN THIS CASE WE USE COL 2
  ## AS THE "base_col"
  base_col <- 1
  if(sum(v[,base_col]) == 0) base_col <- 2
  ## INTIATE A MATRIX TO STORE OUR ADJUSTED v VALUES
  x <- adjv <- matrix(0,ncol=ncol(v),nrow(v))
  ## x IS THE PROBS DIVIDED BY THE BASE COLUMN, AND LOGGED
  x <- log(v/v[,base_col])
  
  ## NOW ADJUST COPYING PROBS BY LAMBDA
  x[,colindex] <- x[,colindex] + lambda
  ## GET OUT OF LOG SPACE
  adjv <- apply(x,2,function(i) exp(i))
  ## RENORMALISE
  adjv <- adjv / rowSums(adjv)
  ## ESTIMATE THE LOG-LIKELIHOOD
  llik <- -(loglik(v=adjv,data=data,nsamps=nsamps,useall=useall,slike=slike) + dnorm(lambda,0,10,log=TRUE))
  return(llik)
}



print(paste0("estimating likelihoods for chromosome: ", chrom, " in ", pop))
n_haps <- ncol(paintedchrom)/10
n_samps <- 10
n_snps <- nrow(snps)
region_ids <- ancreg_list
region_ids2 <- region_ids[colSums(ind_copy_probs)!=0]
n_regs2 <- length(region_ids2)
self_reg <- region_ids[!region_ids%in%region_ids2]
useallsamps <- TRUE

n_regs <- length(ancreg_list)
ind_copy_probs <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
colnames(ind_copy_probs) <- 1:n_regs
ind_vec <- seq(1,ncol(paintedchrom),n_samps)
for(i in 1:n_haps)
{
  print(paste("estimating probs for hap:", i))
  cp <- matrix(0,nrow=n_samps,ncol=n_regs)
  colnames(cp) <- 1:n_regs
  for(j in 1:n_samps)
  {
    tmp_cp <- table(paintedchromreg[,(ind_vec[i]:(ind_vec[i]+9))[j]])
    tmp_cp <- tmp_cp/sum(tmp_cp)
    cp[j,names(tmp_cp)] <- tmp_cp
  }
  ind_copy_probs[((i-1)*10+1):((i-1)*10+10),] <- cp
}

ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)

v <- c()
for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))


mle <- matrix(0,nrow=n_snps,ncol=3*n_regs2)
cnames <- c()
for(i in region_ids2) cnames <- c(cnames,i,i,i)
colnames(mle) <- paste(cnames,rep(c("likelihood","lambda","P"),n_regs2),sep=".")
for(i in as.numeric(iis))
{
  print(i)
  for(reg_index in 1:length(region_ids))
  {
    reg_id <- region_ids[reg_index]
    if(reg_id != self_reg)    
    #if(reg_id == "East_Africa_Nilo-Saharan")
    {
      data <- paintedchromreg[i,]
      ## test if local copier copies local
      painted <- hapregs[copiercopies[i,paintedchrom[i,]]] == self_reg
      data[painted] <- NA
      lambda <- 0
      opt1 <- optim(lambda,par.loglik,
                    data=data,
                    nsamps=n_samps,useall=useallsamps,
                    slike=samp2use,colindex=reg_index,v=v,
                    method="Nelder-Mead")
      ## COMPUTE NULL
      null_lik <- loglik(v,data,nsamps=n_samps,useall=useallsamps,slike=samp2use)
      
      lrt <- (-2*(null_lik) + (2*(-opt1$value)))
      p <- -log10(pchisq(q=lrt,df=1,lower.tail=F))
      mle[i,grep(reg_id,colnames(mle))]  <- c(-opt1$value,opt1$par,p)
    }
  }
}

## is the likelihood ratio test correct here??
## implement the test and run genome-wide?
## need to run the whocopiedcopies script on chroms 11-22


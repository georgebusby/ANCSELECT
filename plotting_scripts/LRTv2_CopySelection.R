###################################################################
###################################################################
##################
## SCRIPT TO GENERATE LOG LIKELIHOOD RATIO TEST FOR INDIVIDUAL SNPS
## COMPARED TO GENOME(CHROMOSOME)-WIDE PAINTINGS
##################
###################################################################
###################################################################

temp <- commandArgs()
# temp <- c("",
#         "FULAI",
#         2,
#         "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalAllChromsPP.samples.out.gz",
#         "/mnt/kwiat/well/human/george/chromopainter2/analysislists/FULAInolocal.idfile.txt",
#         "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalAllChromsPP.likelihoods.gz")
# 
pop <- temp[2]
mainchrom <- temp[3]
in_file <- temp[4]
id_file <- temp[5]
out_file<- temp[6]

library("bigmemory")

###################################################################
###################################################################
## SOME IMPORTANT VARIABLES
options(scipen=999,digits=20)
## INFORMATION ON SNPS - CHROMOSOME,POSITION,ALLELES ETC
## DIRECTORY WITH FILE WITH SNP INFO IN THEM
#snp_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
snp_dir <- "/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
snp_file_pre <- "AllPops330KChrom"
snp_file_pos <- "phased.legend.gz"
## NUMBER OF SAMPLES TO USE TO GENERATE LIKELIHOODS
n_samps <- 10 # IF 1 THEN JUST USES THE FIRST
useallsamps <- TRUE ## USE ALL SAMPLES?
samp2use <- 1 ## IF useallsamps == FALSE, THEN USE THIS SAMPLE
###################################################################
###################################################################
###################################################################
## INFORMATION ON THE SAMPLES
## POPFILE MUST HAVE TWO COLUMNS: "Ethnic_Group" AND "Region"
## THIS FILE WILL BE READ TO GROUP POPULATIONS INTO REGIONS
#pop_file <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
pop_file <- "/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
popkey <- read.table(pop_file,header=T)
###################################################################
###################################################################
###################################################################
### FUNCTIONS ###
###################################################################
###################################################################
###################################################################
### PROGRAM ###
## 00 LOAD PAINTING SAMPLES
## LOAD SAMPLES: NOTE THESE ARE GENERATED USING THE PY-PROG
## AND SHOULD BE SNPS AS COLUMNS AND HAPS AS COLUMNS
#lines3 <- read.big.matrix(in_file,type="char",sep=" ")
lines3 <- as.big.matrix(read.table(in_file,colClasses="integer"))

## 01 GET POPULATION INFO PLUS INFO ON REGIONS, NUMBERS OF HAPS ETC.
ids <- read.table(id_file)    
pops <- as.character(ids[,2])
regions <- sapply(pops,function(x){as.character(popkey$Region[popkey$Ethnic_Group==x])})
ids <- cbind(ids,regions)
n_regs <- length(unique(regions))
region_ids <- levels(popkey$Region)
n_haps <- sum(ids$V2==pop)*2
## FOR EACH DONOR HAPLOTYPE WE WANT TO KNOW THE IDENTITY OF THE RGION THAT IS
## COMES FROM. THIS ALIGNS THE DONORS IN THE *samples.out FILES WITH THEIR
## REGIONAL IDENTITY
donor_hap_vec <- c()
for(i in 1:nrow(ids)) donor_hap_vec <- c(donor_hap_vec,ids$regions[i],ids$regions[i])
pop_hap_vec <- c()
for(i in 1:nrow(ids)) pop_hap_vec <- c(pop_hap_vec,as.character(ids$V2[i]),as.character(ids$V2[i]))
id_hap_vec <- c()
for(i in 1:nrow(ids)) id_hap_vec <- c(id_hap_vec,as.character(ids$V1[i]),as.character(ids$V1[i]))

## 02 LOAD SNP INFO
snps <- c()
for(chrom in 1:22)
{
    if(chrom < 10) chrom  <- paste0("0",chrom)
    snp_file <- paste(snp_dir,snp_file_pre,chrom,snp_file_pos,sep="")
    snp <- read.table(snp_file,header=F)
    snps <- rbind(snps,snp)
}    
colnames(snps) <- c("chrom","rsid","pos","a0","a1")

## 03 LEAVE-ONE-OUT COPYING PROBS
## COMPUTE GENOME-WIDE COPYING-PROBS BASED ON NUMBER OF 
## SNPS COPIED FROM EACH REGION ACROSS ALL CHROMSOMES EXCEPT
## THE ONE THAT WE'RE ANALYSING
## WE AVERAGE THESE ACROSS ALL SAMPLES ????

chrom <- mainchrom
rows <- snps$chrom!=chrom
#lines4 <- lines3[rows,]
print(paste0("estimating likelihoods for chromosome: ", chrom, " in ", pop))
ind_copy_probs <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
colnames(ind_copy_probs) <- 1:n_regs
ind_vec <- seq(1,ncol(lines3),n_samps)
for(i in 1:n_haps)
{
    print(paste("esimating probs for hap:", i))
    cp <- matrix(0,nrow=n_samps,ncol=n_regs)
    colnames(cp) <- 1:n_regs
    for(j in 1:n_samps)
    {
        tmp_cp <- table(donor_hap_vec[lines3[rows,(ind_vec[i]:(ind_vec[i]+9))[j]]])
        tmp_cp <- tmp_cp/sum(tmp_cp)
        cp[j,names(tmp_cp)] <- tmp_cp
    }
    ind_copy_probs[((i-1)*10+1):((i-1)*10+10),] <- cp
}

## THESE ARE THE MAIN COPYING PROPORTIONS
ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)

## REMOVE REGION THAT IS ARE NOT COPIED FROM A REGION ID VECTOR
region_ids2 <- region_ids[colSums(ind_copy_probs)!=0]
n_regs2 <- length(region_ids2)
self_reg <- region_ids[!region_ids%in%region_ids2]

## 04 GET SNPS ON CURRENT CHROMOSOME AND ESTIMATE NULL LIKELIHOODS
snps2 <- snps[snps$chrom==chrom,]
n_snps <- nrow(snps2)
lines4 <- as.big.matrix(lines3[snps$chrom==chrom,])

## DEFINE WHETHER WE AVERAGE COPY PROBS ACROSS ALL SAMPLES OR JUST USE ONE SAMPLE
v <- c()
if(useallsamps == TRUE) for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))
if(useallsamps == FALSE) v <- ind_copy_probs[seq(samp2use,nrow(ind_copy_probs),by=10),]

## NULL SNP LIKELIHOODS
snp_liks <- c()
for(i in 1:nrow(lines4))
{
    ll <- loglik(v,donor_hap_vec[lines4[i,]],nsamps=n_samps,useall=useallsamps,slike=samp2use)
    snp_liks <- c(snp_liks,ll)
}
    
####### NEW STANDARD, COOKBOOK APPROACH

## 00 DEFINE A MATRIX OF MUS
mus <- matrix(0,nr=nrow(ind_copy_probs)/10,nc=ncol(ind_copy_probs))
cnt <- 1
for(i in seq(1,n_haps*n_samps,n_samps))
{
    mus[cnt,] <- colSums(ind_copy_probs[i:(i+9),])/n_samps
    cnt <- cnt+1
}

## LOGIT OF MUS
mus <- log(mus/( 1-mus))

####################
## LOGISTIC FUNCTIONS
find_beta <- function(beta,sum_y,mu){return(sum_y - sum( 1/(1+exp(-(mu+beta)))))}
find_riged_beta <- function(beta,sum_y,mu,sigma_sq){return(sum_y -0.5*(beta^2)/sigma_sq-sum( 1/(1+exp(-(mu+beta)))))}
sigma_sq<-2

####################
## 01 PROGRAM

donor_reg_vec <- (1:length(region_ids))[region_ids!=self_reg]

out_mat <- snps[snps$chrom==mainchrom,]
for(donor_reg in donor_reg_vec)
{
  mu_i <- mus[,donor_reg] 
  ps <- rep(0,sum(snps$chrom==mainchrom))
  lrs <- rep(0,sum(snps$chrom==mainchrom))
#  ps_riged <- rep(0,sum(snps$chrom==2))
  for(snp_index in 1:sum(snps$chrom==mainchrom))
  {
    sum_yi <- sum(donor_hap_vec[lines4[snp_index,]] == donor_reg) / n_samps
    beta_hat<-uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
#    beta_hat_riged<-uniroot(find_riged_beta,interval=c(-10,10),sum_y=sum_yi,mu=mu_i,sigma_sq=sigma_sq)$root

    LRT<-2*(beta_hat*sum_yi*n_samps+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat)))))
#    LRT_riged<-2*(beta_hat_riged*sum_yi*n_samps-n_samps*0.5*(beta_hat_riged^2)/sigma_sq+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat_riged)))))
    lrs[snp_index] <- LRT
     
    ps[snp_index]<- -pchisq(LRT,df = 1,lower.tail = F,log.p = T)/log(10)
#    ps_riged[snp_index]<- -pchisq(LRT_riged,df = 1,lower.tail = F,log.p = T)/log(10)
  }
  out_mat <- cbind(out_mat,lrs,ps)
  colnames(out_mat)[ncol(out_mat)-1] <- paste(region_ids[donor_reg],"likelihood",sep=".")
  colnames(out_mat)[ncol(out_mat)] <- paste(region_ids[donor_reg],"P",sep=".")
  
}


write.table(out_mat,file=out_file,quote=F,col.names=T,row.names=F)


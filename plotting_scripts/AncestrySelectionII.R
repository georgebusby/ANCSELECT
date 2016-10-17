###################################################################
###################################################################
##################
## SCRIPT TO GENERATE LOG LIKELIHOOD RATIO TEST FOR INDIVIDUAL SNPS
## COMPARED TO GENOME(CHROMOSOME)-WIDE PAINTINGS
##################
###################################################################
###################################################################
library("h5")
#library("bigmemory")
library("data.table")

temp <- c("",
        "FULAI",
        "22",
        "/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5",
        "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/ancestry_selection/FULAInolocalChrom22.ancestryselectionII.gz")

temp <- commandArgs()
  
pop <- temp[2]
mainchrom <- temp[3]
datafile <- h5file(temp[4], mode = 'r')
outfile<- temp[5]


popkey_file <- "kwiat/2/bayes/users/george/popgen/analysis3/popgen/ancestry_selection/MalariaGenAdmixturePopulationOverviewNSAA.txt"
popkey_file <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/ancestry_selection/MalariaGenAdmixturePopulationOverviewNSAA.txt"
snp_dir <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
#snp_dir <- "/kwiat/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
pop_file <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
#pop_file <- "/kwiat/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"

options(scipen=999,digits=20)
n_samps <- 10
ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                            "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                            "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)


###################################################################
###################################################################
### FUNCTIONS ###
hap2sampleindex <- function(hap,nsamps=10){
  ## finds the first sample index for a haplotype
  sample <- (hap*nsamps)-(nsamps-1)
  return(sample)
}

###################################################################
###################################################################
## THE MAIN LIKELIHOOD FUNTION:
## THIS TAKES:
##      v: COPYING PROBS PER HAPLOTYPE; THESE ARE ESTIMATED FROM GENOME-WIDE SAMPLES
##      data: A VECTOR OF HAPLOTYPES AT A GIVEN SNP
##      nsamps: THE NUMBER OF SAMPLES OF EACH HAPLOTYPE AT THIS SNP
##      useall: SHOLD ALL PAINTED HAPLOTPYES BE USED TO GENERATE LIKELIHOOD?
##      slike: IF useall = FALSE ; THEN USE THIS SAMPLE TO GENERATE LIKELIHOOD
##      NOTE: THE COLUMNS OF V SHOULD CORRESPOND TO THE COPYING PROBS
##          OF COPYING FROM THE REGION WITH THE SAME INDEX
##          THAT IS, COLUMN 1 OF v SHOULD DESCRIBE THE PROB OF COPYING
##          FROM THE REGION IDENTIFIED BY 1 IN data
##      NOTE: THAT DATA SHOULD BE n_haps * n_samps long

loglik <- function(v,data,nsamps=n_samps)
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
  #if(useall == TRUE ) prob <- sum(log(rowSums(prob)/ncol(prob)))
  ## THIS IS GBJB'S PREVIOUS VERSION, WHICH LOOKS LIKE IT WORKS
  ## BUT IS PERHAPS INCORRECT THING TO DO
  #prob <- mean(apply(prob,2,function(x)sum(log(x), na.rm = T)), na.rm = T)
  prob <- sum(log(rowSums(prob)/n_samps), na.rm = T)
  return(prob)
}
## THIS FUNCTION OPTIMISES THE VALUE OF LAMBDA, BASED ON MLE USING OPTIM
##  RECALL THAT v MUST HAVE n_regions COLUMNS, WITH EACH COLUMN DESCRIBING
##  THE COPYING PROBS FROM A GIVEN REGION INDEX. THESE ARE THE INDICES IN data
##      lambda: IS THE VALUE WE WANT TO OPTIMISE (THE PROPORIONAL INCREASE IN 
##          COPYING FROM A GIVEN DONOR REGION)
##      colindex: THE COLUMN OF v THAT WE WANT TO ADJUST TO GET OPTIMUM VALUE
##          OF LAMBDA - IE THE REGION OF INTEREST

par.loglik <- function(v,data,nsamps=n_samps,lambda,colindex)
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
  llik <- -(loglik(v=adjv,data=data,nsamps=nsamps) + dnorm(lambda,0,10,log=TRUE))
  return(llik)
}
###################################################################
###################################################################
###################################################################
### PROGRAM ###
## ALL DATA IS IN BIG HDF5 FILE
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- readDataSet(datafile[paste0("/paintings/samples/individuals")])
colnames(psamples) <- c("ind","region","X")

####################################################################################
## GET PAINTINGS ACROSS ALL POP OF INTEREST
analysis <- "nonlocal"
psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
# 2 haps per sample!!
psampleshap <- hap2sampleindex(psamplesind,2)
psampleshap <- sort(c(psampleshap,psampleshap+1))
psamplesindsamp <- hap2sampleindex(psampleshap)
tmp <- c()
for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))

## LOAD ALL PAINTINGS
chroms <- 1:22
paintings <- snps <- c()
for(chrom in chroms)
{
  print(paste0("loading chromosome: ",chrom))
  if(chrom<10) chrom <- paste0("0",chrom)
  tmpchrom <- data.table(datafile[paste0("paintings/chrom",chrom,"/",analysis)][,tmp])
  paintings <- rbind(paintings,tmpchrom)
  tmpmap <- data.frame(readDataSet(datafile[paste0("paintings/chrom",chrom,"/map")]),stringsAsFactors = F)
  colnames(tmpmap) <- c("position","recrate")
  tmpmap <- data.frame(apply(tmpmap,2,as.numeric))
  tmpsnps <- data.frame(readDataSet(datafile[paste0("paintings/chrom",chrom,"/snps")]))
  colnames(tmpsnps) <- c("chrom","rsid","pos","a0","a1")
  tmpsnps <- cbind(tmpsnps,tmpmap$recrate)
  snps <- rbind(snps,tmpsnps)  
}

snps <- data.table(snps)
colnames(snps)[ncol(snps)] <- "recrate"
snps$pos <- as.numeric(as.character(snps$pos))
snps$recrate <- as.numeric(as.character(snps$recrate))
n_snps <- sum(snps$chrom==as.numeric(mainchrom))
###################################################################
## 01 GET POPULATION INFO PLUS INFO ON REGIONS, NUMBERS OF HAPS ETC.
paintedchrom <- paintings[snps$chrom==as.numeric(mainchrom),]
n_regs <- length(regions)
n_haps <- ncol(paintedchrom)/n_samps
## SWITCH DONORS TO REGIONS
happops <- c()
for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
hapregs <- c()
for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
hapregs <- factor(hapregs,levels=ancreg_list)
paintedchromreg <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchromreg[,i] <- hapregs[paintedchrom[[i]]]
paintedchrompop <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchrompop[,i] <- happops[paintedchrom[[i]]]
###################################################################
## 03 LEAVE-ONE-OUT COPYING PROBS
## COMPUTE GENOME-WIDE COPYING-PROBS BASED ON NUMBER OF 
## SNPS COPIED FROM EACH REGION ACROSS ALL CHROMSOMES EXCEPT
## THE ONE THAT WE'RE ANALYSING
## THESE ARE THE MAIN COPYING PROPORTIONS
print(paste0("estimating likelihoods for chromosome: ", mainchrom, " in ", pop))
ind_copy_probs <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
colnames(ind_copy_probs) <- regions
rows <- which(snps$chrom!=as.numeric(mainchrom))
dt <- paintings[rows,]
for(j in 1:ncol(dt)) set(dt,j=j,value=hapregs[dt[[j]]])
for(j in 1:ncol(dt))
{
  print(paste0("generating genome-wide copying probs for haplotype/sample: ",j,"/",n_haps*n_samps));
  cptab <- table(dt[[j]])
  ind_copy_probs[j,names(cptab)] <- cptab
}
ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)

## REMOVE REGION THAT IS ARE NOT COPIED FROM A REGION ID VECTOR
region_ids2 <- region_ids[colSums(ind_copy_probs)!=0]
n_regs2 <- length(region_ids2)
self_reg <- region_ids[!region_ids%in%region_ids2]
self_reg2 <- as.character(popkey$Region[popkey$Ethnic_Group==pop])

###################################################################
## 04 GET SNPS ON CURRENT CHROMOSOME AND ESTIMATE NULL LIKELIHOODS
## AVERAGE COPY PROBS ACROSS 10 SAMPLES
v <- c()
for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))

####################################################################
## 05a ORIGINAL GB/CC APPROACH
## GET ID OF WHO COPIER COPIES
copiercopies <- datafile[paste0("/lengths/chrom",mainchrom,"/copiercopies")]

## NOW FIND MLE OF LAMBDA
mle <- matrix(0,nrow=n_snps,ncol=3*n_regs2)
cnames <- c()
for(i in region_ids2) cnames <- c(cnames,i,i,i)
colnames(mle) <- paste(cnames,rep(c("likelihood","lambda","P"),n_regs2),sep=".")
for(i in 1:n_snps)
{
  print(i)
  #pcdone <- signif((i/n_snps)*100,2)
  #if(pcdone%%10 == 0) print(paste(pcdone," % through snps"))
  for(reg_index in 1:length(region_ids))
  {
    reg_id <- region_ids[reg_index]
    if(!reg_id %in% self_reg)
    {
      data <- paintedchromreg[i,]
      ## test if local copier copies local
      painted <- hapregs[copiercopies[i,unlist(paintedchrom[i])]] == self_reg
      data[painted] <- NA
      lambda <- 0
      #opt1 <- optim(lambda,par.loglik,data=data,nsamps=n_samps,colindex=reg_index,v=v,method="Nelder-Mead")
      opt2 <- optimise(par.loglik,interval = c(-5,5),data=data,nsamps=n_samps,colindex=reg_index,v=v)
      
      ## COMPUTE NULL
      null_lik <- -(loglik(v,data,nsamps=n_samps) + dnorm(lambda,0,10,log=TRUE))
      test_lik <- opt2$objective
      test_lam <- opt2$minimum
      lrt <- (2*null_lik) - (2*test_lik)
      p <- -log10(pchisq(q=lrt,df=1,lower.tail=F))
      mle[i,grep(reg_id,colnames(mle))]  <- c(test_lik,test_lam,p)
    }
  }
}

## I WANT TO TRY TO SPEED THIS UP A BIT
## 

# ## DOES COPIER COPY SELF?
# painted <- as.matrix(paintedchrom[1:10])
# p2 <- cbind(sort(rep(1:nrow(painted),ncol(painted))),as.vector(painted))
# p2[,3] <- copiercopies[p2[,1],p2[,2]]
# for(i in 1:ncol(painted))
# {
#   hapregs[copiercopies[1,painted[1,]]]
# }
# copiercopies[1,]
mle <- cbind(snps[snps$chrom==as.numeric(mainchrom),],mle)
##mle <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalChrom22PP.likelihoods.gz",header = T)

all_out <- mle

####################################################################
## 05b RYAN'S NEW STANDARD, COOKBOOK APPROACH
## 00 DEFINE A MATRIX OF MUS
#mus <- matrix(0,nr=nrow(ind_copy_probs)/10,nc=ncol(ind_copy_probs))
# cnt <- 1
# for(i in seq(1,n_haps*n_samps,n_samps))
# {
#     mus[cnt,] <- colSums(ind_copy_probs[i:(i+9),])/n_samps
#     cnt <- cnt+1
# }
## LOGIT OF MUS
mus <- ind_copy_probs
mus <- log(mus/( 1-mus))

####################################################################
## RYAN'S LOGISTIC FUNCTIONS
find_beta <- function(beta,sum_y,mu)
{
  return(sum_y - sum( 1/(1+exp(-(mu+beta)))))
}
find_ridged_beta <- function(beta,sum_y,mu,sigma_sq)
{
  return(sum_y -0.5*(beta^2)/sigma_sq-sum( 1/(1+exp(-(mu+beta)))))
}
####################################################################
sigma_sq <- 2
donor_reg_vec <- (1:length(region_ids))[region_ids!=self_reg]
out_mat <- snps[snps$chrom==as.numeric(mainchrom),]
for(reg_index in 1:length(region_ids))
{
  reg_id <- region_ids[reg_index]
  if(!reg_id %in% self_reg)
  {
    print(paste0("generating likelihoods for: ",reg_id))
    mu_i <-  mus[,reg_index]
    ps <- lrs <- betas <- rep(0,n_snps)
    # ps_riged <- rep(0,sum(snps$chrom==2))
    for(snp_index in 1:n_snps)
    {
      ## NO RECIPROCAL VERSION
      yi <- paintedchromreg[snp_index,] == reg_index
      sum_yi <- sum(yi) / n_samps
      
      ## RECIPROCAL VERSION
      data <- paintedchromreg[snp_index,]
      ## test if local copier copies local
      print("here")
      painted <- hapregs[copiercopies[snp_index,unlist(paintedchrom[i])]] == self_reg
      print("or here")
      data[painted] <- NA
      ## make a mu table for each sampled painting??
      sum_yi <- sum(data == reg_index,na.rm=T)
      mu_i <-  mus[,reg_index][!painted]
      if(sum_yi > 0)
      {
        beta_hat<-uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
        #beta_hat_riged<-uniroot(find_riged_beta,interval=c(-10,10),sum_y=sum_yi,mu=mu_i,sigma_sq=sigma_sq)$root
        LRT<-2*(beta_hat*sum_yi*n_samps+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat)))))
        #LRT_riged<-2*(beta_hat_riged*sum_yi*n_samps-n_samps*0.5*(beta_hat_riged^2)/sigma_sq+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat_riged)))))
        lrs[snp_index] <- LRT
        ps[snp_index]<- -pchisq(LRT,df = 1,lower.tail = F,log.p = T)/log(10)
        betas[snp_index] <- beta_hat
        #ps_riged[snp_index]<- -pchisq(LRT_riged,df = 1,lower.tail = F,log.p = T)/log(10)
      } else
      {
        lrs[snp_index] <- ps[snp_index] <- betas[snp_index] <- NA
      }
    }
    out_mat <- cbind(out_mat,lrs,ps,betas)
    colnames(out_mat)[ncol(out_mat)-2] <- paste(reg_id,"likelihood",sep=".")
    colnames(out_mat)[ncol(out_mat)-1] <- paste(reg_id,"P",sep=".")
    colnames(out_mat)[ncol(out_mat)] <- paste(reg_id,"beta",sep=".")
  }
}

out_mat <- data.table(out_mat)
#old_out <- out_mat

## GENERATE SOME EMPIRICAL P-VALUES
## eg
## tmp <- ecdf("ALL non-target chromosome SNPS -log10 p vlaues")("chromosome -log10 p values"))
## THIS WILL GIVE YOU A LIST OF TOP X REGIONS ETC

# 
# 
# #write.table(out_mat,file=out_file,quote=F,col.names=T,row.names=F)
# 
colnames(out_mat) <- gsub("likelihood","likelihoodII",colnames(out_mat))
colnames(out_mat) <- gsub(".P",".PII",colnames(out_mat))
mle_out <- cbind(mle,out_mat[6:ncol(out_mat)])

####################################################################
## 06 CHRIS'S MVN METHOD
#Set things up
# n_ind <- n_haps / n_samps;
#haps <- matrix(donor_hap_vec[lines4[1:nrow(lines4),]],byrow = F,nr=n_snps)
#Get individuals averages
## AVERAGE ACROSS THE 10 SAMPLES?
# avg <- ind_copy_probs
# avg <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
# colnames(avg) <- 1:n_regs
# ind_vec <- seq(1,ncol(lines4),n_samps)
# for(i in 1:n_haps)
# {
#   print(paste("estimating genome-wide avg for hap:", i))
#   cp <- matrix(0,nrow=n_samps,ncol=n_regs)
#   colnames(cp) <- 1:n_regs
#   for(j in 1:n_samps)
#   {
#     tmp_cp <- table(donor_hap_vec[lines3[rows,(ind_vec[i]:(ind_vec[i]+9))[j]]])
#     tmp_cp <- tmp_cp/sum(tmp_cp)
#     cp[j,names(tmp_cp)] <- tmp_cp
#   }
#   avg[((i-1)*10+1):((i-1)*10+10),] <- cp
# }
# 
# avg <- avg/rowSums(avg)
# 
# #Residual
# res <- paintedchromreg;
# for(i in 1:ncol(res))
# {
#   for(j in 1:length(region_ids))
#   {
#     res[paintedchromreg[,i] == j,i] <- (1 - avg[i,j])
#   }
#   print(i);
# }
# 
# #Get the total residual deviation
# first <- seq(1,n_haps*10,by=n_samps);
# deviant <- array(0,c(n_samps,n_snps,length(region_ids)));
# for(j in 1:n_samps)
# {
#   index <- first + (j-1);
#   for(i in 1:length(region_ids))
#   {
#     tmp <- res[,index];
#     tmp[res[,index] != i] = 0;
#     deviant[j,,i] = rowSums(tmp);
#   }
#   print(j);
# }
# 
# random <- sample(1:n_snps,500);
# tmp <- deviant[,,which(colSums(avg) != 0)];
# x <- array(0,dim(tmp)[-1])
# for(i in 1:n_samps) x <- x + tmp[i,,]
# x <- x / n_samps;
# mu <- colSums(x[random,]) / nrow(x[random,]);
# sigma <- cov(x[random,]);
# # anc = 6;
# # plot(-log10(pnorm(x[,anc],mu[anc],sd=sqrt(sigma[anc,anc]))))
# prior.sigma = array(0,dim(sigma))
# diag(prior.sigma) = 0.0001;
# sigma <- sigma + prior.sigma;
# output <- array(NA,c(n_snps,2));
# for(i in 1:n_snps)
#   output[i,1] <- t(x[i,] - mu) %*% solve(sigma) %*% (x[i,] - mu);
# output[,2] <- pchisq(output[,1],ncol(sigma),lower=FALSE)
# 
# marg.pval = array(NA,dim(x));
# for(i in 1:ncol(marg.pval))
#   marg.pval[,i] = pchisq((x[,i]-mu[i])^2/diag(sigma)[i],1,lower=FALSE);
# 
# colnames(output) = c("MVNchisq","MVNp");
# colnames(x) <- region_ids[!region_ids%in%self_reg]
# colnames(marg.pval) <- paste(region_ids[!region_ids%in%self_reg],".MVNp",sep="")
# res <- data.frame(snps[1:n_snps,],output,x,marg.pval);
# 
# all_out <- cbind(mle_out,res[,6:ncol(res)])

options(digits = 5);

all_out <- mle_out
write.csv(all_out,outfile)







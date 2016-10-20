###################################################################
###################################################################
##################
## SCRIPT TO GENERATE LOG LIKELIHOOD RATIO TEST FOR INDIVIDUAL SNPS
## COMPARED TO GENOME(CHROMOSOME)-WIDE PAINTINGS
##################
###################################################################
###################################################################

temp <- c("",
        "FULAI",
        22,
        "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalAllChromsPP.samples.out.gz",
        "/mnt/kwiat/well/human/george/chromopainter2/analysislists/FULAInolocal.idfile.txt",
        "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalAllChromsPP.likelihoods.gz")

# temp <- commandArgs()
  
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
snp_dir <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
#snp_dir <- "/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
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
pop_file <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
#pop_file <- "/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
popkey <- read.table(pop_file,header=T)
###################################################################
###################################################################
###################################################################
###################################################################
### FUNCTIONS ###
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
  if(useall == TRUE) prob <- mean(apply(prob,2,function(x)sum(log(x))))
  return(prob)
}

## THIS FUNCTION OPTIMISES THE VALUE OF LAMBDA, BASED ON MLE USING OPTIM
##  RECALL THAT v MUST HAVE n_regions COLUMNS, WITH EACH COLUMN DESCRIBING
##  THE COPYING PROBS FROM A GIVEN REGION INDEX. THESE ARE THE INDICES IN data
##      lambda: IS THE VALUE WE WANT TO OPTIMISE (THE PROPORIONAL INCREASE IN 
##          COPYING FROM A GIVEN DONOR REGION)
##      colindex: THE COLUMN OF v THAT WE WANT TO ADJUST TO GET OPTIMUM VALUE
##          OF LAMBDA - IE THE REGION OF INTEREST

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
###################################################################
###################################################################
###################################################################
### PROGRAM ###
###################################################################
## 00 LOAD PAINTING SAMPLES
## LOAD SAMPLES: NOTE THESE ARE GENERATED USING THE PY-PROG
## AND SHOULD BE SNPS AS COLUMNS AND HAPS AS COLUMNS
#lines3 <- read.big.matrix(in_file,type="char",sep=" ")
lines3 <- as.big.matrix(read.table(in_file,colClasses="integer"))

###################################################################
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

###################################################################
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

###################################################################
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

###################################################################
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
  
####################################################################
## 05a ORIGINAL GB/CC APPROACH
## NOW FIND MLE OF LAMBDA
mle <- matrix(0,nrow=n_snps,ncol=3*n_regs2)
cnames <- c()
for(i in region_ids2) cnames <- c(cnames,i,i,i)
colnames(mle) <- paste(cnames,rep(c("likelihood","lambda","P"),n_regs2),sep=".")
for(i in 1:n_snps)
{
  pcdone <- signif((i/n_snps)*100,2)
  if(pcdone%%10 == 0) print(paste(pcdone," % through snps"))
  for(reg_index in 1:length(region_ids))
  {
    reg_id <- region_ids[reg_index]
    if(reg_id != self_reg)
    {
      lambda <- 0
      opt1 <- optim(lambda,par.loglik,
                    data=donor_hap_vec[lines4[i,]],
                    nsamps=n_samps,useall=useallsamps,
                    slike=samp2use,colindex=reg_index,v=v,
                    method="Nelder-Mead")
      lrt <- (-2*(snp_liks[i]) + (2*(-opt1$value)))
      p <- -log10(pchisq(q=lrt,df=1,lower.tail=F))
      mle[i,grep(reg_id,colnames(mle))]  <- c(-opt1$value,opt1$par,p)
    }
  }
}


mle <- cbind(snps2,snp_liks,mle)

#mle <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalChrom22PP.likelihoods.gz",header = T)

####################################################################
## 05b RYAN'S NEW STANDARD, COOKBOOK APPROACH
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
## RYAN'S LOGISTIC FUNCTIONS
find_beta <- function(beta,sum_y,mu){
  return(sum_y - sum( 1/(1+exp(-(mu+beta)))))
  }

find_ridged_beta <- function(beta,sum_y,mu,sigma_sq){
  return(sum_y -0.5*(beta^2)/sigma_sq-sum( 1/(1+exp(-(mu+beta)))))
}

sigma_sq <- 2

donor_reg_vec <- (1:length(region_ids))[region_ids!=self_reg]
out_mat <- snps[snps$chrom==mainchrom,]
for(reg_index in 1:length(region_ids))
{
  reg_id <- region_ids[reg_index]
  if(reg_id != self_reg)
  {
    print(paste0("generating likelihoods for: ",reg_id))
    mu_i <- mus[,reg_index] 
    ps <- lrs <- rep(0,n_snps)
    # ps_riged <- rep(0,sum(snps$chrom==2))
    for(snp_index in 1:n_snps)
    {
      sum_yi <- sum(donor_hap_vec[lines4[snp_index,]] == reg_index) / n_samps
      if(sum_yi > 0)
      {
        beta_hat<-uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
        #beta_hat_riged<-uniroot(find_riged_beta,interval=c(-10,10),sum_y=sum_yi,mu=mu_i,sigma_sq=sigma_sq)$root
    
        LRT<-2*(beta_hat*sum_yi*n_samps+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat)))))
        #LRT_riged<-2*(beta_hat_riged*sum_yi*n_samps-n_samps*0.5*(beta_hat_riged^2)/sigma_sq+n_samps*sum(log((1+exp(mu_i))/(1+exp(mu_i+beta_hat_riged)))))
        lrs[snp_index] <- LRT
         
        ps[snp_index]<- -pchisq(LRT,df = 1,lower.tail = F,log.p = T)/log(10)
        #ps_riged[snp_index]<- -pchisq(LRT_riged,df = 1,lower.tail = F,log.p = T)/log(10)
      } else
      {
        lrs[snp_index] <- ps[snp_index] <- NA
      }
    }
    out_mat <- cbind(out_mat,lrs,ps)
    colnames(out_mat)[ncol(out_mat)-1] <- paste(reg_id,"likelihood",sep=".")
    colnames(out_mat)[ncol(out_mat)] <- paste(reg_id,"P",sep=".")
  }
}


#write.table(out_mat,file=out_file,quote=F,col.names=T,row.names=F)

colnames(out_mat) <- gsub("likelihood","likelihoodII",colnames(out_mat))
colnames(out_mat) <- gsub(".P",".PII",colnames(out_mat))

mle_out <- cbind(mle,out_mat[,6:ncol(out_mat)])

####################################################################
## 06 CHRIS'S MVN METHOD
#Set things up
n_ind <- n_haps / n_samps;
haps <- matrix(donor_hap_vec[lines4[1:nrow(lines4),]],byrow = F,nr=n_snps)
#Get individuals averages

## AVERAGE ACROSS THE 10 SAMPLES?
avg <- array(NA, c(n_haps,length(region_ids)));
for(i in 1:n_haps)
{
  k <- 1+((i-1)*10)
  index <- k:(k+9)
  for(j in 1:length(region_ids)) avg[i,j] <- sum(haps[,index] == j)/n_samps;  
}

avg <-  avg / n_snps;

#Residual
res <- haps;
for(i in 1:n_haps)
{
  k <- 1+((i-1)*10)
  index <- k:(k+9)
  for(j in 1:length(region_ids))
  {
    res[res[,index] == j] <- (1 - avg[i,j])
  }
  print(i);
}

#Get the total residual deviation
first <- seq(1,n_haps*10,by=n_samps);
deviant <- array(0,c(n_samps,n_snps,length(region_ids)));
for(j in 1:n_samps)
{
  index = first + (j-1);
  for(i in 1:length(region_ids))
  {
    tmp <- res[,index];
    tmp[haps[index,] != i] = 0;
    deviant[j,,i] = rowSums(tmp);
  }
  print(j);
}

#x11(type="Xlib",width=15);
#barplot(t(counts[1,,]),beside=FALSE,border=NA,space=0,col=2:7);
random <- sample(1:n_snps,500);
tmp <- deviant[,,which(colSums(avg) != 0)];
x <- array(0,dim(tmp)[-1])
for(i in 1:n_samps) x <- x + tmp[i,,] 
x <- x / n_samps;
mu <- colSums(x[random,]) / nrow(x[random,]);
sigma <- cov(x[random,]);
# anc = 6;
# plot(-log10(pnorm(x[,anc],mu[anc],sd=sqrt(sigma[anc,anc]))))
prior.sigma = array(0,dim(sigma))
diag(prior.sigma) = 0.0001;
sigma <- sigma + prior.sigma;
output <- array(NA,c(n_snps,2));
for(i in 1:n_snps)
  output[i,1] <- t(x[i,] - mu) %*% solve(sigma) %*% (x[i,] - mu);
output[,2] <- pchisq(output[,1],ncol(sigma),lower=FALSE)

marg.pval = array(NA,dim(x));
for(i in 1:ncol(marg.pval))
  marg.pval[,i] = pchisq((x[,i]-mu[i])^2/diag(sigma)[i],1,lower=FALSE);

colnames(output) = c("MVNchisq","MVNp");
colnames(x) <- region_ids[region_ids!=self_reg]
colnames(marg.pval) <- paste(region_ids[region_ids!=self_reg],".MVNp",sep="")
res <- data.frame(snps[1:n_snps,],output,x,marg.pval);

all_out <- cbind(mle_out,res[,6:ncol(res)])
options(digits = 5);


write.csv(all_out,out_file);







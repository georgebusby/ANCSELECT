###################################################################
###################################################################
##################
## SCRIPT TO GENERATE LOG LIKELIHOOD RATIO TEST FOR INDIVIDUAL SNPS
## COMPARED TO GENOME(CHROMOSOME)-WIDE PAINTINGS
##################
###################################################################
###################################################################
library(rhdf5)
library(data.table)

# temp <- c("",
#         "AFAR",
#         "22",
#         "/mnt/kwiat/data/1/galton/users/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5")

temp <- commandArgs()
pop <- temp[2]
mainchrom <- temp[3]
#datafile <- H5Fopen(temp[4])

paste0 <- function(...) {
  paste(...,sep="")
}


# main_dir <- "/kwiat/1/galton/users/george/copy_selection/"
# main_dir <- "/mnt/kwiat/data/1/galton/users/george/copy_selection/"
main_dir <- "/well/malariagen/malariagen/human/george/copy_selection2/copy_selection/"
# main_dir <- "/mnt/kwiat/well/human/george/copy_selection2/copy_selection/"

outfile <- paste0(main_dir,"output/",pop,"nolocalChrom",mainchrom,".ancestryselectionV.gz")

# popkey_file <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/ancestry_selection/MalariaGenAdmixturePopulationOverviewNSAA.txt"
# snp_dir <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
# pop_file <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"

popkey_file <- paste0(main_dir,"analysis_lists/MalariaGenAdmixturePopulationOverviewNSAA.txt")
snp_dir <- paste0(main_dir,"snpfiles/")
pop_file <- paste0(main_dir,"analysis_lists/populationOverviewCopyProbs.txt")
hdf5_file <- paste0(main_dir,"hdf5files/MalariaGenSelectionPaintings.hdf5")
datafile <- H5Fopen(hdf5_file)


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


###################################################################
###################################################################
###################################################################
### PROGRAM ###
## ALL DATA IS IN BIG HDF5 FILE
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- t(h5read(datafile,paste0("/paintings/samples/individuals")))
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
snps <- c()
for(chrom in chroms)
{
  print(paste0("loading chromosome: ",chrom))
  if(chrom<10) chrom <- paste0("0",chrom)
  tmpchrom <- data.table(t(h5read(datafile,paste0("paintings/chrom",chrom,"/",analysis),index = list(tmp,NULL))))
  if(chrom == "01")
  {
    paintings <- tmpchrom
  } else
  {
    paintings <- rbind(paintings,tmpchrom)
  }
  
  tmpmap <- data.frame(t(h5read(datafile,paste0("paintings/chrom",chrom,"/map"))),stringsAsFactors = F)
  colnames(tmpmap) <- c("position","recrate")
  tmpmap <- data.frame(apply(tmpmap,2,as.numeric))
  tmpsnps <- data.frame(t(h5read(datafile,paste0("paintings/chrom",chrom,"/snps"))))
  colnames(tmpsnps) <- c("chrom","rsid","pos","a0","a1")
  tmpsnps <- cbind(tmpsnps,tmpmap$recrate)
  snps <- rbind(snps,tmpsnps)  
}
rm(tmpchrom)

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
rm(dt)
rm(paintings)
ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)

## REMOVE REGION THAT IS ARE NOT COPIED FROM A REGION ID VECTOR
region_ids2 <- regions[colSums(ind_copy_probs)!=0]
n_regs2 <- length(region_ids2)
self_reg <- regions[!regions%in%region_ids2]
self_reg2 <- as.character(popkey$Region[popkey$Ethnic_Group==pop])

###################################################################
## 04 GET SNPS ON CURRENT CHROMOSOME AND ESTIMATE NULL LIKELIHOODS
## AVERAGE COPY PROBS ACROSS 10 SAMPLES
v <- c()
for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))

####################################################################
## 05a ORIGINAL GB/CC APPROACH
## GET REGION OF WHO COPIER COPIES
copiercopies <- data.table(t(h5read(datafile,paste0("/lengths/chrom",mainchrom,"/copiercopies"))))
for(j in 1:ncol(copiercopies)) set(copiercopies,j=j,value=hapregs[copiercopies[[j]]])
paintedchromreg <- data.table(paintedchromreg)

## SET EVERYTHING UP
mle <- matrix(0,nrow=n_snps,ncol=(n_regs2*7)+1)
cnames <- c()
for(i in region_ids2) cnames <- c(cnames,rep(i,7))
colnames(mle) <- c("pc.drop",paste(cnames,
                                   rep(c("prop","GB.lik","GB.beta","GB.P",
                                         "RC.lik","RC.beta","RC.P"),n_regs2),sep="."))

## LOGIT OF MUS
mus <- ind_copy_probs
mus <- log(mus/( 1-mus))

### STORE A DATAFRAME FOR THE MVN TEST
pchrom <- data.frame(paintedchromreg)

for(i in 1:n_snps)
{
  print(paste("generating likelihoods for snp:",i, "/", n_snps))
  #pcdone <- signif((i/n_snps)*100,2)
  #if(pcdone%%10 == 0) print(paste(pcdone," % through snps"))
  data <- unlist(paintedchromreg[i])
  ## test if local copier copies local
  copied_haps <- paste0("V",paintedchrom[i])
  copied_hapsreg <- as.character(hapregs[unlist(paintedchrom[i])])
  #painted <- copiercopies[i=i,j=copied_haps, with = F] == self_reg
  painted <- copiercopies[i=i,j=copied_haps, with = F] != copied_hapsreg
  data[painted] <- NA
  pchrom[i,] <- data
  perc.dropped <- sum(painted)/length(painted)
  mle[i,1] <- perc.dropped
  new_props <- table(data)
  prop_cols <- paste0(ancreg_list[as.numeric(names(new_props))],".prop")
  mle[i,prop_cols] <- table(data)
  for(reg_index in 1:length(regions))
  {
    
    reg_id <- regions[reg_index]
    if(!reg_id %in% self_reg)
    {
      ###################################################
      ## GEORGE'S HACKED LRT TEST 
      lambda <- 0
      #opt1 <- optim(lambda,par.loglik,data=data,nsamps=n_samps,colindex=reg_index,v=v,method="Nelder-Mead")
      opt2 <- optimise(par.loglik,interval = c(-5,5),data=data,nsamps=n_samps,colindex=reg_index,v=v)
      ## COMPUTE NULL
      null_lik <- -(loglik(v,data,nsamps=n_samps) + dnorm(lambda,0,10,log=TRUE))
      test_lik <- opt2$objective
      beta <- opt2$minimum
      lrt <- (2*null_lik) - (2*test_lik)
      p <- -log10(pchisq(q=lrt,df=1,lower.tail=F))
      gb_cols <- paste0(reg_id,c(".GB.lik",".GB.beta",".GB.P"))
      mle[i,gb_cols]  <- c(lrt,beta,p)
      ###################################################
      ## RYAN'S PRINCIPALLED VERSION
      sum_yi <- sum(data == reg_index,na.rm=T)
      mu_i <-  mus[,reg_index][!painted]
      avg_num_samps <- n_samps/(nrow(ind_copy_probs)/sum(!is.na(data)))
      if(sum_yi > 0)
      {
        BETA <- uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
        sum_logmu <- sum(log((1+exp(mu_i))/(1+exp(mu_i+BETA))))
        beta_y <- BETA*sum_yi
        LRT <- 2*(beta_y + sum_logmu)/avg_num_samps ## this might be controversial
        P <- -log10(pchisq(q=LRT,df=1,lower.tail=F))
      } else
      {
        BETA <- LRT <- P <- NA
      }
      rc_cols <- paste0(reg_id,c(".RC.lik",".RC.beta",".RC.P"))
      mle[i,rc_cols] <- c(LRT,BETA,P)
    }
  }
}

mle <- cbind(snps[snps$chrom==as.numeric(mainchrom),],mle)
all_out <- data.table(mle)

## GENERATE SOME EMPIRICAL P-VALUES
## eg
## tmp <- ecdf("ALL non-target chromosome SNPS -log10 p vlaues")("chromosome -log10 p values"))
## THIS WILL GIVE YOU A LIST OF TOP X REGIONS ETC
####################################################################
## 06 CHRIS'S MVN METHOD
## Set things up
n_ind <- n_haps / n_samps;
avg <- ind_copy_probs
res <- pchrom
## CALCULATE RESIDUALS
res[is.na(res)] <- 0
for(i in 1:ncol(res))
{
  for(j in 1:length(regions))
  {
    res[res[,i] == j,i] <- (1 - avg[i,j])
  }
  print(i);
}

## GET THE TOTAL RESIDUAL DEVIATION 
first <- seq(1,n_haps*10,by=n_samps);
deviant <- array(0,c(n_samps,n_snps,length(regions)));
for(j in 1:n_samps)
{
  index <- first + (j-1);
  for(i in 1:length(regions))
  {
    tmp <- res[,index];
    tmp[pchrom[,index] != i] = 0;
    deviant[j,,i] = rowSums(tmp);
  }
  print(j);
}

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
output[,2] <- -log10(pchisq(output[,1],ncol(sigma),lower=FALSE))

marg.pval = array(NA,dim(x));
for(i in 1:ncol(marg.pval))
  marg.pval[,i] = -log10(pchisq((x[,i]-mu[i])^2/diag(sigma)[i],1,lower=FALSE));

colnames(output) <- c("MVNchisq","MVNp");
colnames(x) <- paste(regions[!regions%in%self_reg],".MVNprops",sep = "")
colnames(marg.pval) <- paste(regions[!regions%in%self_reg],".MVNp",sep="")
output <- data.frame(output,x,marg.pval);

all_out <- cbind(all_out,output)
options(digits = 5);
write.table(all_out, file = outfile, row.names = F, quote = F, sep = ",")

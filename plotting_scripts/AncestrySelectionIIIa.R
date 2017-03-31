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
#        "AFAR",
#        "22",
#        "/mnt/kwiat/data/1/galton/users/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5")

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
### PROGRAM ###
## ALL DATA IS IN BIG HDF5 FILE
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- t(h5read(datafile,paste0("/paintings/samples/individuals")))
colnames(psamples) <- c("ind","region","X")

pop <- "FULAI"
#for(pop in pops[2:5])
#{
  ####################################################################################
  psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
  # 2 haps per sample!!
  psampleshap <- hap2sampleindex(psamplesind,2)
  psampleshap <- sort(c(psampleshap,psampleshap+1))
  psamplesindsamp <- hap2sampleindex(psampleshap)
  tmp <- c()
  for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
  n_haps <- length(tmp)
  n_regs <- length(regions)
  happops <- c()
  for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
  happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
  hapregs <- c()
  for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
  hapregs <- factor(hapregs,levels=ancreg_list)
  analysis <- "nonlocal"
  
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
    
    # #########################################
    # ## test if local copier copies local
    # copiercopies <- t(h5read(datafile,paste0("/lengths/chrom",chrom,"/copiercopies")))
    # copiercopies <- data.table(copiercopies)
    # for(j in 1:ncol(copiercopies)) set(copiercopies,j=j,value=hapregs[copiercopies[[j]]]==hapregs[j])
    # copymat <- as.matrix(copiercopies, byrow = T)
    # indices <- unlist(tmpchrom)
    # indices <- cbind(rep(1:nrow(tmpchrom),n_haps),indices)
    # indices <- cbind(indices,copymat[cbind(as.numeric(indices[,1]),as.numeric(indices[,2]))])
    # copymat <- matrix(indices[,3],nr=nrow(copymat),nc = ncol(tmpchrom), byrow = F)
    # ## copymat tells us whether each hap should be masked at that SNP
    # if(chrom == "01")
    # {
    #   maskings <- copymat
    # } else
    # {
    #   maskings <- rbind(maskings,copymat)
    # }
    tmpmap <- data.frame(t(h5read(datafile,paste0("paintings/chrom",chrom,"/map"))),stringsAsFactors = F)
    colnames(tmpmap) <- c("position","recrate")
    tmpmap <- data.frame(apply(tmpmap,2,as.numeric))
    tmpsnps <- data.frame(t(h5read(datafile,paste0("paintings/chrom",chrom,"/snps"))))
    colnames(tmpsnps) <- c("chrom","rsid","pos","a0","a1")
    tmpsnps <- cbind(tmpsnps,tmpmap$recrate)
    snps <- rbind(snps,tmpsnps)  
  }
  
  ###################################################################
  ## 02 GET SNP INFO ETC.
  ## GET PAINTINGS ACROSS ALL POP OF INTEREST
  snps <- data.table(snps)
  colnames(snps)[ncol(snps)] <- "recrate"
  snps$pos <- as.numeric(as.character(snps$pos))
  snps$recrate <- as.numeric(as.character(snps$recrate))
  
  for(mainchrom in 1:22)
  {
    print(paste0("estimating likelihoods for chromosome: ", mainchrom, " in ", pop))
    if(mainchrom < 10) mainchrom <- paste0("0",mainchrom)
    
    n_snps <- sum(snps$chrom==as.numeric(mainchrom))
    paintedchrom <- paintings[snps$chrom==as.numeric(mainchrom),]
    n_haps <- ncol(paintedchrom)/n_samps
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
    ind_copy_probs <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
    colnames(ind_copy_probs) <- regions
    rows <- which(snps$chrom!=as.numeric(mainchrom))
    dt <- data.frame(paintings[rows,])
    #dt2 <- data.frame(maskings[rows,])
    for(j in 1:ncol(dt)) set(dt,j=j,value=hapregs[dt[[j]]])
    for(j in 1:ncol(dt))
    {
      #print(paste0("generating genome-wide copying probs for haplotype/sample: ",j,"/",n_haps*n_samps));
      cptab <- table(dt[[j]])
      #cptab <- table(unlist(dt[[j]])[unlist(dt2[[j]])==0])
      ind_copy_probs[j,names(cptab)] <- cptab
    }
    
    ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)
    pchrom <- paintedchromreg
    ## REMOVE REGION THAT IS ARE NOT COPIED FROM A REGION ID VECTOR
    region_ids2 <- regions[colSums(ind_copy_probs)!=0]
    n_regs2 <- length(region_ids2)
    self_reg <- regions[!regions%in%region_ids2]
    self_reg2 <- as.character(popkey$Region[popkey$Ethnic_Group==pop])
    
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
     # print(j);
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
    
    all_out <- cbind(snps[snps$chrom==as.numeric(mainchrom),],output)
    options(digits = 5);
    outfile <- paste0(main_dir,"output/",pop,"nolocalChrom",mainchrom,".ancestryselectionIIIa.gz")
    write.table(all_out, file = outfile, row.names = F, quote = F, sep = ",")
  }
#}



# rm(dt)
# rm(paintings)


###################################################################
# ## 04 GET SNPS ON CURRENT CHROMOSOME AND ESTIMATE NULL LIKELIHOODS
# ## AVERAGE COPY PROBS ACROSS 10 SAMPLES
# v <- c()
# for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))

####################################################################
## 05a ORIGINAL GB/CC APPROACH
## GET REGION OF WHO COPIER COPIES
## WE DON'T NEED THIS ANYMORE AS WE'VE STORED THIS INFO ABOVE
# copiercopies <- data.table(t(h5read(datafile,paste0("/lengths/chrom",mainchrom,"/copiercopies"))))
# for(j in 1:ncol(copiercopies)) set(copiercopies,j=j,value=hapregs[copiercopies[[j]]])


## SET EVERYTHING UP
# mle <- matrix(0,nrow=n_snps,ncol=(n_regs2*7)+1)
# cnames <- c()
# for(i in region_ids2) cnames <- c(cnames,rep(i,7))
# colnames(mle) <- c("pc.drop",paste(cnames,
#                                    rep(c("prop","GB.lik","GB.beta","GB.P",
#                                          "RC.lik","RC.beta","RC.P"),n_regs2),sep="."))
# ## LOGIT OF MUS
# mus <- ind_copy_probs
# mus <- log(mus/( 1-mus))
# 
# ## STORE A DATAFRAME FOR THE MVN TEST
# pchrom <- paintedchromreg
# ## TURN INTO DATATABLE FOR SPEED
# paintedchromreg <- data.table(paintedchromreg)
# 
# for(i in 1:n_snps)
# {
#   print(paste("generating likelihoods for snp:",i, "/", n_snps))
#   data <- unlist(paintedchromreg[i,])
#   ## test if local copier copies local: NOW DONE ABOVE
#   #copied_haps <- paste0("V",paintedchrom[i])
#   #copied_hapsreg <- as.character(hapregs[unlist(paintedchrom[i])])
#   #painted <- copiercopies[i=i,j=copied_haps, with = F] == self_reg
#   #painted <- copiercopies[i=i,j=copied_haps, with = F] != copied_hapsreg
#   #painted <- maskings[i,] == 0
#   #data[painted] <- NA
#   
#   pchrom[i,] <- data
#   #perc.dropped <- sum(painted)/length(painted)
#   perc.dropped <- 0
#   mle[i,1] <- perc.dropped
#   new_props <- table(data)
#   prop_cols <- paste0(ancreg_list[as.numeric(names(new_props))],".prop")
#   mle[i,prop_cols] <- table(data)
#   for(reg_index in 1:length(regions))
#   {
#     reg_id <- regions[reg_index]
#     if(!reg_id %in% self_reg)
#     {
#       ###################################################
#       ## GEORGE'S HACKED LRT TEST 
#       lambda <- 0
#       #opt1 <- optim(lambda,par.loglik,data=data,nsamps=n_samps,colindex=reg_index,v=v,method="Nelder-Mead")
#       opt2 <- optimise(par.loglik,interval = c(-10,10),data=data,nsamps=n_samps,colindex=reg_index,v=v)
#       ## COMPUTE NULL
#       null_lik <- -(loglik(v,data,nsamps=n_samps) + dnorm(lambda,0,10,log=TRUE))
#       test_lik <- opt2$objective
#       beta <- opt2$minimum
#       lrt <- (2*null_lik) - (2*test_lik)
#       p <- -log10(pchisq(q=lrt,df=1,lower.tail=F))
#       gb_cols <- paste0(reg_id,c(".GB.lik",".GB.beta",".GB.P"))
#       mle[i,gb_cols]  <- c(lrt,beta,p)
#       ###################################################
#       ## RYAN'S PRINCIPALLED VERSION
#       sum_yi <- sum(data == reg_index,na.rm=T)
#       #mu_i <-  mus[,reg_index][!painted]
#       mu_i <-  mus[,reg_index]
#       avg_num_samps <- n_samps/(nrow(ind_copy_probs)/sum(!is.na(data)))
#       if(sum_yi > 0)
#       {
#         BETA <- uniroot(find_beta,c(-10,10),sum_y=sum_yi,mu=mu_i)$root
#         sum_logmu <- sum(log((1+exp(mu_i))/(1+exp(mu_i+BETA))))
#         beta_y <- BETA*sum_yi
#         LRT <- 2*(beta_y + sum_logmu)/avg_num_samps ## this might be controversial
#         P <- -log10(pchisq(q=LRT,df=1,lower.tail=F))
#       } else
#       {
#         BETA <- LRT <- P <- NA
#       }
#       rc_cols <- paste0(reg_id,c(".RC.lik",".RC.beta",".RC.P"))
#       mle[i,rc_cols] <- c(LRT,BETA,P)
#     }
#   }
# }
# 
# mle <- cbind(snps[snps$chrom==as.numeric(mainchrom),],mle)
# all_out <- data.table(mle)

## GENERATE SOME EMPIRICAL P-VALUES
## eg
## tmp <- ecdf("ALL non-target chromosome SNPS -log10 p vlaues")("chromosome -log10 p values"))
## THIS WILL GIVE YOU A LIST OF TOP X REGIONS ETC

###
### SCRIPT TO INFER INDIVIDUAL ANCESTRY PROPORTIONS FOR ALL INDS ###

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"ANCSELECT/"))
###########################################################
## DEFINE DATAFILES
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"

## LOAD POPKEY FILE ##
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                            "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                            "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

analysis <- "nonlocal"

##################################################################################
## LOAD INDIVIDUAL CLUSTER ASSIGNMENT:: NEEDS SERVER ACCESS
popplot <- popkey$Ethnic_Group
final_clusts <- vector("list",length(popplot))
for(i in 1:length(popplot))
{
  ii <- as.character(popplot[i])
  names(final_clusts)[i] <- ii
  if(ii =="SEMI.BANTU") ii <- "SEMI-BANTU"
  iinds <- scan(paste0("/mnt/kwiat/home/popgen/scripts/finalpoplists/",ii,"finalinds.txt"),what="char")
  final_clusts[[i]] <- iinds
}
paintedinds <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/samplelists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP.inds",sep=" ")

## CHANGE ACTUAL SAMPLE IDS TO CP LABELS
final_clusts2<- lapply(final_clusts,function(x){as.character(paintedinds[match(unlist(x),paintedinds[,2]),1])})
##################################################################################

## nnls functions
library(nnls)
source("~/repos/ANCSELECT/plotting_scripts/globeNNLS.R")

## directories
in_dir <- "/mnt/kwiat/well/human/george/chromopainter2/output/"
length_dir <- "/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/globetrotter/input/"


alllocrecmat <- allnonrecmat <- c()
for(i in 1:60)
{
  pop <- popkey$Ethnic_Group[i]
  reg <- popkey$RegionM[i]
  if(reg == "East_Africa_Afroasiatic" ) reg <- "East_Africa_Nilo-Saharan"
  print(paste0("getting individual NNLS results for individuals from ",pop))
  lengths_in <- paste0(length_dir,"MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinal",reg,".chunklengths.out")
  lengths <- read.table(lengths_in,header = T, row.names = 1, as.is = T)
  colnames(lengths) <- gsub("\\.","\\-",colnames(lengths))
  if(reg == "Eurasia")    lengths <- t(rowsAsMapClusts(final_clusts2,t(lengths)))
  ## MAKE DONOR POP COPYING VECTORS
  predmat <- rowsAsMapClusts(final_clusts2,lengths,mean)
  predmat <- predmat/rowSums(predmat)
  colnames(predmat) <- gsub("\\.","\\-",colnames(predmat))
  
  inds <- unique(unlist(final_clusts2[pop]))
  localrecmat <- nonlocalrecmat <- matrix(0,nc = length(popplot),nr = length(inds))
  colnames(localrecmat) <- colnames(nonlocalrecmat) <- popplot
  rownames(localrecmat) <- rownames(nonlocalrecmat) <- inds  
  
  for(ind in inds)
  {
    rec <- lengths[ind,colnames(predmat)]
    rec <- rec/sum(rec)
    ## INFER FOR ALL POTENTIAL DONORS
    matpred <- predmat[rownames(predmat)!=pop,]
    mod <- getoverallfit(matpred,rec)$x
    mod[mod<0.01] <- 0
    mod <- mod/sum(mod)
    localrecmat[ind,names(mod)] <- mod
    
    ## NOW FOR ALL NON-LOCAL SDONORS
    localpops <- popkey$Ethnic_Group[popkey$RegionM==reg]
    
    matpred <- predmat[!rownames(predmat)%in%localpops,]
    mod <- getoverallfit(matpred,rec)$x
    mod[mod<0.01] <- 0
    mod <- mod/sum(mod)
    nonlocalrecmat[ind,names(mod)] <- mod
  }
  
  alllocrecmat <- rbind(alllocrecmat,localrecmat)
  allnonrecmat <- rbind(allnonrecmat,nonlocalrecmat)

}


write.table(alllocrecmat,file = "data/MalariaGenAdmixtureLocalAncestryProps.txt",
            col.names = T, quote = F, row.names = T)

write.table(allnonrecmat,file = "data/MalariaGenAdmixtureNonLocalAncestryProps.txt",
            col.names = T, quote = F, row.names = T)

##################################################
## FOR JACK-KNIFING; SOME DERIVATIVE OF THIS, I
## HAVEN'T RUN THROUGH THE CODE YET ...
# chroms <- 1:22
# for(chrom in chroms)
# {
#   print(chrom)
#   if(chrom < 10) chrom <- paste0("0",chrom)
#   tmp <- read.table(paste0(in_dir,pop,"nolocalChrom",chrom,".chunklengths.out"),
#                     header = T, as.is=T, row.names = 1)
#   tmp <- 
#     if(chrom == "01")
#     {
#       lengths <- tmp
#     } else 
#     {
#       lengths <- lengths + tmp
#     }
# }
##################################################


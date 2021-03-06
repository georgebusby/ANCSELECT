#######################################################
### SCRIPT TO WORK WITH COPYING PROBS FOR GWAS INDS ###
#######################################################

## THE FILES ARE STORED AS HDF5 FILES ON /well AND THIS
## SCRIPT HELPS TO PULL OUT THE RELEVANT INFO

library(h5) ## INSTALL IF NECESSARY

#######################################################
## WE'LL WORK ON A SINGLE CHROMOSOME, COPYPROBS ARE SPLIT ACROSS
## POULATIONS WITH 22 FILES, ONE PER CHROMOSOME
## DEFINE SOME VARIABLES ##
chrom <- "09"
pop <- "Kenya"
in_dir <- "/mnt/kwiat/well/human/george/GWAS_painting/"

samplefile <- paste0(in_dir,"data/Superset2_GWAS-2.5M_b37.PCA.sample")
ethnicfile <- paste0(in_dir,"data/CP1_clin_phenotypes_ALL_25JAN2016.csv")
popsamplefile <- paste0(in_dir,"abacus/",pop,"/ALL/",pop,"GWASphasedChrom",chrom,".sample")


## in_dir IS THE RELEVANT DIRECTORY MOUNTED ON MY OWN COMPUTER -- THIS FULL PATH IS BELOW
#in_dir <- "/well/malariagen/malariagen/human/george/GWAS_painting/chromopainter/hdf5/"
hdf5file <- paste0(in_dir,"chromopainter/hdf5/",pop,"Chrom",chrom,"Copyprobsperlocus.h5")
haplofile <- paste0(in_dir,"chromopainter/hdf5/",pop,"Chrom",chrom,"PhasedHaplotypes.h5")

#######################################################
## DEFINE OUR DATAFILES
datafile <- h5file(hdf5file, mode = "r")
hapsfile <- h5file(haplofile, mode = "r")
#######################################################
## THE BASIC LAYOUT OF datafile IS:
##    /copyprobsperlocus/chrom04/
##    /copyprobsperlocus/chrom04/donors  = painting donors
##    /copyprobsperlocus/chrom04/map  = info on snps
## THEN WE HAVE THE ACTUAL DATA STORED IN CHUNKS OF 200 HAPS
## TO GET AT THE INDS, WE LOOK HERE
##    /copyprobsperlocus/chrom04/haps/1-100 etc
## THE ACTUAL COPYING PROBS ARE HERE
##    /copyprobsperlocus/chrom04/probs/1-100 etc
## LET'S TRY WITH A SINGLE CHROMOSOME: 22
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
# 1 /copyprobsperlocus/chrom09/haps/1-100       1              1 HAP 1 5423088011_R01C01
# 2 /copyprobsperlocus/chrom09/haps/1-100       2              2 HAP 2 5423088011_R01C01
# 3 /copyprobsperlocus/chrom09/haps/1-100       3              3 HAP 1 5423088011_R02C01
# 4 /copyprobsperlocus/chrom09/haps/1-100       4              4 HAP 2 5423088011_R02C01
# 5 /copyprobsperlocus/chrom09/haps/1-100       5              5 HAP 1 5423088011_R03C01
# 6 /copyprobsperlocus/chrom09/haps/1-100       6              6 HAP 2 5423088011_R03C01

## WHICH TELLS ME THAT HAPLOTYPE 1 IS HAP1 OF IND 5423088011_R01C01 AND IS STORED A POSITION 1 IN hdfdataset
##                               2 IS HAP2 OF IND 5423088011_R01C01 AT POSITION 2, etc
## LET'S GET THE CASE/CONTROL STATUS
samples <- read.table(samplefile, header = T)
samples <- samples[-1,]
hap_casecontrol <- samples$caseorcontrol[match(sapply(as.character(haps$haplotypeID),
                                                      function(x){strsplit(x, split=" ")[[1]][3]}),
                                               samples$gtc_id)]
hap_scl <- samples$sc_sl[match(sapply(as.character(haps$haplotypeID),
                                                      function(x){strsplit(x, split=" ")[[1]][3]}),
                                               samples$gtc_id)]
hap_scl <- as.character(hap_scl)
## LET'S GET ETHNICITY
ethnics <- read.csv(ethnicfile,header = T)
hap_ethnicity <- as.character(ethnics$curated_ethnicity[match(hap_scl,ethnics$sc_sl)])
hap_ethnicity[is.na(hap_ethnicity)] <- "OTHER"
## TO GET INDICES OF FULA, FOR EXAMPLE DO THIS:
# ethnic <- "FULA"
# ethnic_indices <- which(hap_ethnicity==ethnic)
## COMBINE HAPS, CASECONTROL AND ETHNICITY INTO ONE BIG FILE
haps <- cbind(haps,hap_casecontrol,hap_ethnicity)
## THESE FILES CONTAIN ALL INDIVIDUALS IN THE GWAS; 
## GAVIN GENERATED AN UNRELATED SUBSET WHICH CAN BE IDENTIFIED
## BECAUSE THEY WEREN'T USED TO GENERATE PCs
## WE'RE GOING TO USE THE SAMPLES THAT WE HAVE VITERBI
## PAINTINGS FOR
haps_keep <- read.table(popsamplefile, header = T, stringsAsFactors = F, as.is = T)
haps_keep <- haps_keep[-1,]
## THIS MAKES A LIST OF INDIVIDUALS WHO WE WANT TO PAINT WITH
## ABACUS VITERBI METHOD
haps2keep <- sapply(as.character(haps$haplotypeID),
                    function(x){strsplit(x,split = " ")[[1]][3]})
haps2keep <- haps2keep%in%haps_keep$ID_1
#######################################################
## WE WANT TO MAKE SURE THAT WE JUST USE THESE HAPLOTYPES ##
## FOR FUTURE ANALYSIS

#######################################################
## LET'S ALSO FIND THE ORIGINAL HAPLOTYPES
## THESE ARE STORED IN haplofile
## EG FOR THE FIRST 100 INDS IN THE COPYPROBS FILE (200 HAPS)
hapsdataset <- as.character(haps$hdf5dataset[1])
hapsdataset <- gsub("copyprobsperlocus","phased",hapsdataset)  
hapsdataset <- gsub("haps","haplotypes",hapsdataset)
tmp_haps <- data.frame(readDataSet(hapsfile[hapsdataset]))

#######################################################
## 01 LOOK AT SOME SNPS ACROSS CASES AND CONTROLS
snp_index <- 1:nrow(snps)
## THE PROBS ARE STORED AS A 3D ARRAY: n_haps * n_snps * n_dons
probs_available <- list.datasets(datafile[paste0("/copyprobsperlocus/chrom",chrom,"/probs/")])

## DONORS ARE SPLIT INTO POPULATIONS;
## WE WANT TO SUM ACROSS POPS FROM SAME ANCESTRY REGION
## WE CAN GET THE DONORS HERE
donors <- datafile[paste0("/copyprobsperlocus/chrom",chrom,"/donors")][]
## WE CAN CONVERT THE DONORS TO ANCESTRY REGIONS LIKE THIS
popkey <- read.table("data/PopulationKey.txt", header = T, as.is = T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
donor_regions <- popkey$AncestryRegion[match(donors,popkey$Ethnic_Group)]
regions <- sort(unique(popkey$AncestryRegion))

#######################################################
## 02 SO, LET'S GET THE PROBS AT ALL SNPS FOR 
## ALL GAMBIAN CASES
selection <- haps[hap_casecontrol=="CASE",]
probs_haps <- matrix(0,nr=0,nc=length(donors)) ## STORE THE PROBABILITIES OT EACH HAPLOTYPE
probs_snps <- matrix(0,nr=nrow(snps),nc=length(donors)) ## STORE THE CHROMOSOME-WIDE PROBS AT EACH SNP
for(i in as.character(unique(selection$hdf5dataset)))
{
  hap_index <- selection$hdf5dataset_ID[selection$hdf5dataset==i]
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

## AND THERE YOU GO ...



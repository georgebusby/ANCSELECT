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
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
datafile <- h5file('/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5')
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
#####
## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
  ## finds the first sample index for a haplotype
  sample <- (hap*nsamps)-(nsamps-1)
  return(sample)
}


hlapops <- c(popkey$Ethnic_Group[popkey$RegionM=="Western_Africa_Niger-Congo"][8])
chrom <- "06"
copy_reg <- "East_Africa_Nilo-Saharan"

## GET MAP AND POSITION INFO
map <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/map")]),stringsAsFactors = F)
colnames(map) <- c("position","recrate")
map <- data.frame(apply(map,2,as.numeric))
snps <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/snps")]))
colnames(snps) <- c("chrom","rsid","pos","a0","a1")
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- readDataSet(datafile[paste0("/paintings/samples/individuals")])
colnames(psamples) <- c("ind","region","X")
chromlength <- as.numeric(as.character(map$recrate))
chrompos <- as.numeric(as.character(map$position))
chromposI <- c(diff(chrompos),0)

## MAKE A TABLE OF THE HAPLOTYPE WHO ARE BEING COPIED MORE AT THIS REGION
hlacopy <- c()
maxcopy <- c()

pdf("figures/WollofFigure3HLA.pdf", height = 5, width = 10)
#x11()
plot_matrix <- matrix(c(2,5,
                        2,5,
                        2,5,
                        1,5,
                        1,4,
                        1,4,
                        7,4,
                        7,4,
                        7,6,
                        3,6,
                        3,6,
                        3,6), 12,2,byrow = T)
layout(plot_matrix,
       heights = c(rep(1,12)))

plot_range <- c(30e6,35e6)

pop <- "WOLLOF"
print(pop)
############################################
## PLOT LRT TEST
pops <- read.table('/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP.idfile.txt')
pops <- as.character(levels(pops[,2]))

chrom2 <- as.numeric(chrom)
in_dir <- '/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/'
in_file <- paste0(in_dir,pop,'nolocalChrom',chrom,'.ancestryselectionIII.gz')
in_file2 <- paste0(in_dir,pop,'nolocalChrom',chrom,'.ancestryselectionIIIa.gz')
liks <- read.csv(in_file,header=T,as.is=T)
pcols <- grep("RC.P",colnames(liks))

# in_dir <- '/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/'
# in_file <- paste0(in_dir,pop,'nolocalChrom',chrom2,'PP.likelihoods.gz')
# liks <- read.table(in_file,header=T,as.is=T)
# pcols <- grep(".P",colnames(liks))
plot_region <- liks$pos>=plot_range[1]&liks$pos<=plot_range[2]
liks <- liks[plot_region,]
ymax <- max(liks[,pcols])    
ylims <- c(0, 20)    
## FIND THE RANGE OVER WHICH THE PVAL IS HIGH
maxp <- which.max(apply(liks[,pcols],2,max))
maxp <- c(which.max(liks[,pcols[maxp]]), maxp)
pvalue <- liks[maxp[1],pcols[maxp[2]]]
minp <- 3

i <- pvalue
pos1 <- maxp[1]
while(i > minp)
{
  pos1 <- pos1 - 1
  i <- liks[pos1,pcols[maxp[2]]]
}

i <- pvalue
pos2 <- maxp[1]
while(i > minp)
{
  pos2 <- pos2 + 1
  i <- liks[pos2,pcols[maxp[2]]]
}

gene_region1 <- c(liks$pos[pos1],liks$pos[pos2])

par(mar=c(3,6,1,1))
plot(0,0,
     xlim=range(plot_range), type = "n", ylim = ylims,
     axes = F, ylab = expression(-log[10]~italic(P)), xlab = "",
     xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)

for(i in pcols)
{
  reg <- gsub("\\.","-",gsub("\\.RC.P","",colnames(liks)[i]))
  linecol <- pcolshex[regions == reg]
  points(liks$pos,liks[,i], col = linecol, type = "S", lwd = 2)
}
yat <- pretty(ylims)
axis(2, at = yat, labels = yat, las = 2, xpd = T)
xat <- pretty(liks$pos)
axis(1,at=xat,labels = xat/1e6, xpd = T)
mtext(2,at=max(ylims),text="b",adj=0, cex = 1.5, las = 2, line = 5)  

####################################################################################
## SHOW PAINTINGS ACROSS ALL POP OF INTEREST
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
## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
plot_region <- map$position>=plot_range[1]&map$position<=plot_range[2]
paintedchromreg <- paintedchromreg[plot_region,]
paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
for(i in 1:nrow(paintedchromreg))
{
  #print(i)
  tmp <- table(paintedchromreg[i,])
  paintedchromregprop[i,as.numeric(names(tmp))] <- tmp
}
  
chromlength <- as.numeric(as.character(map$recrate))
chrompos <- as.numeric(as.character(map$position))
chromposI <- c(diff(chrompos),0)
chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
chromcols <- pcolshex[chromcolsbreaks]
chromplot <- paintedchromreg

####################################################################################
## PROPORTIONS    
chrompos <- map$position[plot_region]
chromposI <- c(diff(chrompos),0)
chromplot <- paintedchromregprop[,8:1]
colnames(chromplot) <- rev(ancreg_list)
chromplot <- chromplot/rowSums(chromplot)

############################################
### PLOT BARPLOT
par(mar=c(0,6,0.5,1))
barplot(t(chromplot),
        width=chromposI,
        col=rev(pcolshex),xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        xlim=c(0,sum(chromposI)),horiz=F,
        xlab="", ylab = paste0(pop,"\nancestry"))
## add a plot on top to add lines
par(new = T)
plot(0,0,axes = F, ylim=c(0,1), type = "n",
     xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)

mtext(2,at=1,text="a",adj=0, cex = 1.5, las = 2, line = 5)  
############################################
# #######################################################
### PLOT GENES
source("~/repos/glycophorins/external_software/plot_genes.R")
genes = load.genes( "/mnt/kwiat/data/1/galton/malariagen/human/reference/genome-mysql.cse.ucsc.edu/2015-08-18/UCSC_hg19_2015-08-18_refGene.tsv" )
genes <- genes[genes$cdsStartStat!="unk",]
for(plot_run in 1)
{
  if(plot_run == 1) gene_region <- gene_region1
  if(plot_run == 2) gene_region <- gene_region2
  
  par(mar=c(4,6,1,1))
  plot.genes( chromosome = chrom, region = range(gene_region), genes, xaxt = "n" , plot.ylab = "Genes")
  x_at <- pretty(gene_region)
  axis(1,at=x_at,labels=x_at/1e6, xpd = F)
  mtext(1,text = paste("position on chromosome",as.numeric(chrom)), line = 3)
  
  ## PLOT CONNECTING LINES
  region <- range(pretty(gene_region))
  w = which( genes$chromosome == chrom & genes$txEnd >= region[1] & genes$txStart <= region[2] )
  local.genes = genes[w,]
  local.genes$y = NA ;
  local.genes$y[1] = 1 ;
  if( nrow( local.genes ) > 1 ) {
    spacer = ( region[2] - region[1] ) / 10 ;
    maxes = ( local.genes[1,]$txEnd + spacer )
    for( i in 2:nrow( local.genes )) {
      for( l in 1:length( maxes )) {
        if( local.genes$txStart[i] >= maxes[l] ) {
          local.genes$y[i] = l ;
          maxes[l] = local.genes$txEnd[i] + spacer ;
          break ;
        }
      }
      if( is.na( local.genes$y[i] )) {
        maxes = c( maxes, local.genes$txEnd[i] + spacer )
        local.genes$y[i] = length( maxes ) ;
      }
    }
  }
  ymax <- max( 2, max(local.genes$y)+0.5 )
  
  ## add a plot on top 
  par(new = T)
  plot(0,0,axes = F, ylim=c(0,ymax), type = "n",
       xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
  
  ##################################################
  ## ADD LINES
  ## function to find the position of the edge of the device
  line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(0:1, 'inches', 'user'))
    y_off <- diff(grconvertY(0:1, 'inches', 'user'))
    switch(side,
           `1` = par('usr')[3] - line * y_off * lh,
           `2` = par('usr')[1] - line * x_off * lh,
           `3` = par('usr')[4] + line * y_off * lh,
           `4` = par('usr')[2] + line * x_off * lh,
           stop("side must be 1, 2, 3, or 4", call.=FALSE))
  }
  
  segments(plot_range[1], ymax, x1 = gene_region[1], y1 = line2user(1,3), col = "black" ,xpd = T)
  segments(plot_range[2], ymax, x1 = gene_region[2], y1 = line2user(1,3), col = "black" ,xpd = T)
}
mtext(2,at=ymax,text="d",adj=0, cex = 1.5, las = 2, line = 5)  

##########################################################################
## PLOT NUMBER OF DONORS OVER PROPORTION OF ANCESTRY
par(mar=c(3,6,1,1))
plot(0,0,
     xlim=range(plot_range), type = "n", ylim = c(0,1),
     axes = F, ylab = "", xlab = "",
     xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)
abline(h=0.1,col = "grey10", lty = 2)
abline(h=0.05,col = "grey10", lty = 3)
axis(2,labels = NA, tck = 0)
axis(2,at = c(0.25,0.5,0.75,1),las=2)
axis(2,at = 0.1, labels = "0.10", las = 2, lty = 2)
#axis(2,at = 0.05, las = 2, lty = 3)

for(reg in ancreg_list)
{
  doncol <- paste0(gsub("\\-","\\.",reg),".num.dons")
  propcol <- paste0(gsub("\\-","\\.",reg),".prop")
  if(doncol%in%colnames(liks))
  {
    linelwd <- 1
    if(reg == "East_Africa_Nilo-Saharan") linelwd <- 2
    points(liks$pos,liks[,doncol]/liks[,propcol],
           col = pcolshex[ancreg_list==reg], type = "S", lwd = linelwd)
  }
}
xat <- pretty(liks$pos)
axis(1,at=xat,labels = xat/1e6, xpd = T)
mtext(2,at=0.75,text="f",adj=0, cex = 1.5, las = 2, line = 5)  
mtext(3,
      text = "Proportion of unique donors from each donor ancestry across region",
      cex = 0.8)
#########################################################################
## GET DATA FOR WHOLE GENOME
in_root <- paste0("/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/",pop,"nolocalChrom")
in_suff <- ".ancestryselectionIII.gz"

if(!exists(pl))
{  
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
}
##################################################################################
## PLOT HISTOGRAM OF GENOME-WIDE PROP OF UNIQUE DONORS
par(mar=c(3,6,2,1))
gwprops <- pl$East_Africa_Nilo.Saharan.num.dons/pl$East_Africa_Nilo.Saharan.prop

hist(gwprops,
     breaks = 50, main = "", axes = F,
     ylab = "",xlab ="",#  "Genome-wide distribution\nof the proportion of\nunique Nilo-Saharan donors",
     xlim=c(0,1), xaxs = "i", yaxs = "i", col = pcolshex[ancreg_list == "East_Africa_Nilo-Saharan"])
axis(1, at = c(0,0.25,0.5,0.75,1))
axis(2,las = 2)
mtext(2,at=35000,text="e",adj=0, cex = 1.5, las = 2, line = 5)  

## OVERLAY THIS REGION TO SHOW IT'S AN OUTLIER
hist_reg <- liks$pos>gene_region[1]&liks$pos <= gene_region[2]
hist_reg <- median(liks$East_Africa_Nilo.Saharan.num.dons[hist_reg]/liks$East_Africa_Nilo.Saharan.prop[hist_reg])
abline(v=hist_reg, col = "red")
gene_region_text <- round(gene_region1/1e6,1)
empp <- sum(gwprops <= hist_reg)/length(gwprops)
legend("topleft", legend = paste0("Average prop.\nacross ", gene_region_text[1],"-",gene_region_text[2],"Mb\n(red line)\nP =", round(empp,5)), bty = "n")
mtext(3,
      text = "Genome-wide distribution of the prop. of unique Nilo-Saharan donors",
      cex = 0.8, adj = 0.5)
#legend("bottomright",legend = "Genome-wide distribution of the\nproportion unique Nilo-Saharan donors",bty= "n")
#####################################################################
## GET GENOME-WIDE COPYING PROPS FOR THE ANUAK
alllocrecmat <- read.table("data/MalariaGenAdmixtureLocalAncestryProps.txt", 
                           header = T, row.names = 1, sep = " ")
allnonrecmat <- read.table("data/MalariaGenAdmixtureNonLocalAncestryProps.txt", 
                           header = T, row.names = 1)
colnames(allnonrecmat) <- colnames(alllocrecmat) <- gsub("\\.","\\-",colnames(alllocrecmat))

## MAKE A LIST OF CLUSTERS AND THE POPS WITHIN
ancreg_clusters <- vector("list",length(ancreg_list))
for(i in 1:length(ancreg_clusters))
{
  names(ancreg_clusters)[i] <- ancreg_list[i]
  ancreg_clusters[[i]] <- popkey$Ethnic_Group[popkey$RegionM == ancreg_list[i]]
}

allnonrecmat <- t(rowsAsMapClusts(ancreg_clusters,t(allnonrecmat),sum))
alllocrecmat <- t(rowsAsMapClusts(ancreg_clusters,t(alllocrecmat),sum))
  
####################################################################  
####################################################################


####################################################################################
## SHOW PAINTINGS ACROSS ALL POP OF INTEREST
pop <- "WOLLOF"
chrom <- "06"
analysis <- "nonlocal"

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

#####
## SWITCH DONORS TO INDS ...
hapinds <- c()
for(i in 1:nrow(psamples)) hapinds <- c(hapinds,rep(as.character(psamples[i,"ind"]),2))
paintedchromind <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchromind[,i] <- hapinds[paintedchrom[,i]]

reg <- "Western_Africa_Niger-Congo"
paintedchromregional <- paintedchromind
for(i in 1:ncol(paintedchromregional)) paintedchromregional[,i] <- allnonrecmat[paintedchromind[,i],reg]
paintedchromregional <- apply(paintedchromregional,2,as.numeric)

## TO DO : CONVERT THE PAINTED CHROMS TO THE AMOUNT EACH COPIES FROM
## THE REGION OF RELEVANCE AND THEN LOOK TO SEE HOW THE AVERAGE/MAX
## COPYING ALONG THE CHROMOSOME CHANGES
## THIS IS LIKE LOOKING AT COPYIER COPYING, BUT IS SLIGHTLY DIFFERENT
## AS WE'RE LOOKING AT THE ANCESTRY PROPORTIONS OF THE COPYIER
## NB: CHECK THAT WE DON'T HAVE EURASIA NON-LOCAL COPYING PROPS SOMEWHERE ...


dim(paintedchromind)
tmp <- paintedchromind[plot_region,]

whoscopied <- howmanycopied <- matrix(0,nc = 5,nr = nrow(tmp))
for(i in 1:nrow(tmp))
{
  #print(i)
  tmp2 <- sort(table(tmp[i,]), decreasing = T)[1:5]
  whoscopied[i,] <- names(tmp2)
  howmanycopied[i,] <- tmp2
}

###
par(mar=c(3,6,1,1))
plot(liks$pos,howmanycopied[,1],axes = F, type = "S",
     xlab = "",ylab = "",
     xaxs = "i", yaxs = "i", lwd = 2)
mtext(3,text = "Number of times most copied haplotype copied", cex = 0.8)
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)
xat <- pretty(liks$pos)
axis(1,at=xat,labels = xat/1e6, xpd = T)
axis(2,las = 2)

for(i in 1:(length(liks$pos)-1))
{
  xx <- c(liks$pos[i]:liks$pos[(i+1)],liks$pos[(i+1)]:liks$pos[(i)])
  yy <- c(rep(0,length(xx)/2),seq(howmanycopied[i,1],howmanycopied[(i+1),1], length=length(xx)/2))
  polycol <-  popkey$Ethnic_Group==names(final_clusts2)[unlist(lapply(final_clusts2,function(x,y=whoscopied[i,1]){y%in%x}))]
  polycol <- pcolshex[ancreg_list==popkey$RegionM[popkey$Ethnic_Group == names(final_clusts2)[polycol]]]
  polygon(xx,yy,col = polycol, border = NA)
}

mtext(2,at=300,text="g",adj=0, cex = 1.5, las = 2, line = 5)  


#######################################
## NEW PLOT TO SHOW COPIER ANCESTRY TO REPLACE ANUAK 11 PLOT ...
hapinds2 <- c()
#for(i in 1:nrow(psamples)) hapinds2 <- c(hapinds2,paste0(psamples[i,"ind"],":1"),paste0(psamples[i,"ind"],":2"))
for(i in 1:nrow(psamples)) hapinds2 <- c(hapinds2,
                                         hap2sampleindex(which(psamples[,"ind"] == psamples[i,"ind"]),2),
                                         hap2sampleindex(which(psamples[,"ind"] == psamples[i,"ind"]),2)+1)

paintedchromind2 <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
for(i in 1:ncol(paintedchrom)) paintedchromind2[,i] <- hapinds2[paintedchrom[,i]]

## WANT TO GENERATE A NEW BARPLOT WHERE I COMPUTE
## ANCESTRY PROPS AT EACH SNP FOR EACH SEPARATE SET
## OF DONORS THAT ARE BEING COPIED
## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
analysis <- "local"
paintedchrom <- datafile[paste0("/paintings/chrom",chrom,"/",analysis)]

paintedchromprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchrom))
colnames(paintedchromprop) <- ancreg_list
for(i in which(plot_region))
{
  # GET HAPS THAT ARE BEING COPIED
  hapsbeingcopied <- unique(sort(paintedchromind2[i,]))
  # ONLY USE THOSE THAT ARE COPYING THE DEVIANT ANCESTRY?
  hapsbeingcopied<- hapsbeingcopied[which(hapregs[hapsbeingcopied] == copy_reg)]
  
  tmp <- c()
  for(j in hapsbeingcopied) tmp <- c(tmp,(hap2sampleindex(j,10)):(hap2sampleindex(j,10)+9))
  tmp <- paintedchrom[i,tmp]
  tmp <- table(hapregs[tmp])
  paintedchromprop[i,which(names(tmp)%in%colnames(paintedchromprop))] <- tmp
}

paintedchromprop <- paintedchromprop[which(plot_region),]

chromlength <- as.numeric(as.character(map$recrate))
chrompos <- as.numeric(as.character(map$position))
chromposI <- c(diff(chrompos),0)
chromcolsbreaks <- sort(unique(unlist(apply(paintedchromprop,2,unique))))
chromcols <- pcolshex[chromcolsbreaks]
chromplot <- paintedchromprop
####################################################################################
## PROPORTIONS    
## plot order
chromplotorder <- ancreg_list #c(ancreg_list[!ancreg_list%in%reg],reg)
chromplot <- paintedchromprop[,rev(chromplotorder)]
chromplot <- chromplot/rowSums(chromplot)
colnames(chromplot) <- chromplotorder
plotcols <- rev(pcolshex[order(factor(chromplotorder,levels=ancreg_list))])

par(mar=c(0,6,0.5,1))
barplot(t(chromplot),
        width=chromposI[plot_region],
        col=plotcols,xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        horiz=F,#xlim=c(0,sum(chromposI)),
        xlab="", ylab = paste0("Ancestry prop.\nof donor when Nilo-Saharan"),
        xaxs = "i", yaxs = "i")
par(new = T)
plot(0,0,axes = F, ylim=c(0,1), type = "n",
     xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)
mtext(2,at=1,text="c",adj=0, cex = 1.5, las = 2, line = 5)  

###########################################
### PLOT PAINTINGS FOR THE ANUAK IND THAT EVERY COPIES
## HAP 1257
# psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
# psamplesind <- (1:nrow(psamples))[psamples[,"ind"] == "ANUAK11"]
# # 2 haps per sample!!
# psampleshap <- hap2sampleindex(psamplesind,2)
# psampleshap <- sort(c(psampleshap,psampleshap+1))
# psamplesindsamp <- hap2sampleindex(psampleshap)
# tmp <- c()
# for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
# paintedchrom <- datafile[paste0("/paintings/chrom",chrom,"/local")]
# paintedchrom <- paintedchrom[,tmp]
# ## SWITCH DONORS TO REGIONS
# happops <- c()
# for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
# happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
# hapregs <- c()
# for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
# hapregs <- factor(hapregs,levels=ancreg_list)
# paintedchromreg <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
# for(i in 1:ncol(paintedchrom)) paintedchromreg[,i] <- hapregs[paintedchrom[,i]]
# paintedchrompop <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
# for(i in 1:ncol(paintedchrom)) paintedchrompop[,i] <- happops[paintedchrom[,i]]
# 
# ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
# paintedchromreg <- paintedchromreg[plot_region,]
# paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
# for(i in 1:nrow(paintedchromreg))
# {
#   #  print(i)
#   tmp <- table(paintedchromreg[i,])
#   paintedchromregprop[i,as.numeric(names(tmp))] <- tmp
# }
# 
# chromlength <- as.numeric(as.character(map$recrate))
# chrompos <- as.numeric(as.character(map$position))
# chromposI <- c(diff(chrompos),0)
# chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
# chromcols <- pcolshex[chromcolsbreaks]
# chromplot <- paintedchromreg
# ####################################################################################
# ## PROPORTIONS
# chromplot <- paintedchromregprop[,8:1]
# chromplot <- chromplot/rowSums(chromplot)
# colnames(chromplot) <- rev(ancreg_list)

# par(mar=c(0,6,0.5,1))
# barplot(t(chromplot),
#         width=chromposI[plot_region],
#         col=rev(pcolshex),xaxs="i",yaxs="i",
#         space=0,axes=F,xaxt="n",border=NA,
#         horiz=F,#xlim=c(0,sum(chromposI)),
#         xlab="", ylab = paste0("ANUAK 11\nancestry"),
#         xaxs = "i", yaxs = "i")
# 
# # xat <- pretty(c(0,sum(chromposI[plot_region])))
# # xlab <- pretty(chrompos[plot_region])/1e6
# # axis(1, at = xat[1:length(xlab)], labels = xlab)
# ## add a plot on top to add lines
# par(new = T)
# plot(0,0,axes = F, ylim=c(0,1), type = "n",
#      xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
# abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)
# mtext(2,at=1,text="c",adj=0, cex = 1.5, las = 2, line = 5)



dev.off()

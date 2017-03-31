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
chrom <- "04"
copy_reg <- "South_Africa_KhoeSan"

##########################################################
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
##########################################################
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

pdf("figures/Figure4LDB2.pdf", height = 5, width = 10)
plot_matrix <- matrix(c(2,3,
                        1,3,
                        5,4), 3,2,byrow = T)
layout(plot_matrix,
       heights = c(2,2,1),
       widths = c(4,6))

plot_range <- c(14e6,19e6)
hitpops <- c("SUDANESE","SOMALI","ANUAK","GUMUZ","ARI","AFAR","TYGRAY","AMHARA","WOLAYTA","OROMO")
pop <- "SOMALI"
hit1 <- "kgp6561116"
hit2 <- "rs4303974"

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

plot_region <- liks$pos>=plot_range[1]&liks$pos<=plot_range[2]
liks <- liks[plot_region,]
ymax <- max(liks[,pcols])    
ylims <- c(0, 20)    
## FIND THE RANGE OVER WHICH THE PVAL IS HIGH
maxp <- which.max(apply(liks[,pcols],2,max))
maxp <- c(which.max(liks[,pcols[maxp]]), maxp)
pvalue <- liks[maxp[1],pcols[maxp[2]]]
minp <- 0.5

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

gene_region1 <- c(16.4e6,17e6)

par(mar=c(3,6,1,1))
plot(0,0,
     xlim=range(plot_range), type = "n", ylim = ylims,
     axes = F, ylab = expression(-log[10]~italic(P)), xlab = "",
     xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = T)
legend("topright",legend = "Evidence for deviation\nin SOMALI",bty = "n")

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

# ####################################################################################
# ## SHOW PAINTINGS ACROSS ALL POPS OF INTEREST
# psamplesind <- (1:nrow(psamples))[psamples[,"region"] %in% hitpops]
# # 2 haps per sample!!
# psampleshap <- hap2sampleindex(psamplesind,2)
# psampleshap <- sort(c(psampleshap,psampleshap+1))
# psamplesindsamp <- hap2sampleindex(psampleshap)
# tmp <- c()
# for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
# paintedchrom <- datafile[paste0("/paintings/chrom",chrom,"/",analysis)]
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
# ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
# plot_region <- map$position>=plot_range[1]&map$position<=plot_range[2]
# paintedchromreg <- paintedchromreg[plot_region,]
# paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
# for(i in 1:nrow(paintedchromreg))
# {
#   #print(i)
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
# 
####################################################################################
## PROPORTIONS    
chrompos <- map$position[plot_region]
chromposI <- c(diff(chrompos),0)
chromplot <- paintedchromregprop[,8:1]
colnames(chromplot) <- rev(ancreg_list)
chromplot <- chromplot/rowSums(chromplot)

############################################
### PLOT BARPLOT
par(mar=c(2,6,0.5,1))
barplot(t(chromplot),
        width=chromposI,
        col=rev(pcolshex),xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        xlim=c(0,sum(chromposI)),horiz=F,
        xlab="", ylab = paste0("Ancestry across all non\nNiger-Congo speakers from\nEast Africa [n = ",length(psamplesind),"]"))
## add a plot on top to add lines
par(new = T)
plot(0,0,axes = F, ylim=c(0,1), type = "n",
     xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
abline(v=gene_region1, col = "grey10", lty = 1, xpd = F)

abline(v=as.numeric(as.character(snps$pos[snps$rsid==hit1])),col = "red")
axis(1,at=as.numeric(as.character(snps$pos[snps$rsid==hit1])),
     col = "red",labels = hit1)
abline(v=as.numeric(as.character(snps$pos[snps$rsid==hit2])),col = "blue")
axis(1,at=as.numeric(as.character(snps$pos[snps$rsid==hit2])),
     col = "blue",labels = hit2)


mtext(2,at=1,text="a",adj=0, cex = 1.5, las = 2, line = 5)  
############################################
## PLOT PAIRWISE COPYING MATRIX AT SNP OF INTEREST

copyhit4 <- c()
for(hit in c(hit1,hit2))
{
  paintedhit <- paintedchrom[snps$rsid==hit,]
  
  copyhit1 <- matrix(0,nr=length(psamplesind), nc = nrow(psamples)*2)
  rownames(copyhit1) <- psamples[psamplesind,"ind"]
  cnames <- c()
  for(i in psamples[,"ind"]) cnames <- c(cnames,paste0(i,":1"),paste0(i,":2"))
  colnames(copyhit1) <- cnames
  
  for(i in 1:length(psamplesind))
  {
    sampindex <- hap2sampleindex(i,20)
    sampindex <- sampindex:(sampindex+19)
    for(j in paintedhit[sampindex]) copyhit1[i,j] <- copyhit1[i,j] + 1
  }

  ## group all column
  copyhit3 <- copyhit1
  cnamevec <- sapply(colnames(copyhit1),function(x){strsplit(x,split = ":")[[1]][1]})
  tmp <- matrix(0,nr=nrow(copyhit3),nc=ncol(copyhit3)/2)
  rownames(tmp)<- rownames(copyhit3)
  colnames(tmp) <- unique(cnamevec)
  for(i in unique(cnamevec)) tmp[,i] <- rowSums(copyhit3[,cnamevec==i])
  copyhit3 <- tmp
  copyhit3 <- t(rowsAsMapClusts(final_clusts2,t(copyhit3),sum))
  copyhit3 <- rowsAsMapClusts(final_clusts2[hitpops],copyhit3,sum)
  copyhit3 <- copyhit3[rev(rownames(copyhit3)),]
  copyhit4 <- rbind(copyhit4,copyhit3)
}
  
copyhit4 <- copyhit4[,colSums(copyhit4)>0]  
#copyhit4 <- copyhit4[order(factor(rownames(copyhit4),levels = unique(rownames(copyhit4)))),]


  
library(RColorBrewer)
heatcols <- brewer.pal(9,"YlOrRd")
heatcols <- c("white",heatcols)

copyhit4 <- copyhit4/rowSums(copyhit4)
par(mar=c(0.5,2,10,6))
image(1:ncol(copyhit4),
      1:nrow(copyhit4),
      t(copyhit4),
      xlab = "", ylab = "",
      axes = F,
      col = heatcols)
box()
axis(4,at = 1:nrow(copyhit4), label = rownames(copyhit4), lwd = 0, las = 2)
abline(h=10.5,xpd = T)
axis(2,at = 5.25,labels = hit1, col = "red", lwd = 1)
axis(2,at = 15.75,labels = hit2, col = "blue", lwd = 1)

popcol1 <- "X"
for(i in 1:ncol(copyhit4))
{
  pop1 <- colnames(copyhit4)[i]
  popcol <- pcolshex[ancreg_list == popkey$RegionM[popkey$Ethnic_Group==pop1]]
  if(popcol!=popcol1) abline(v=i-0.5, xpd = F, col = "grey10")
  if(pop1 == "GUIGHANAKGAL") pop1 <- "GUIGANA"
  axis(3,at = i,label = tidyNames(pop1,fula = T),
       lwd = 0, las = 2,
       col.axis = popcol)
  popcol1 <- popcol
}

mtext(2,at=nrow(copyhit4)+9.5,text="d",adj=0, cex = 1.5, las = 2, line = 2,xpd = T)  
mtext(3,line = 6.5,cex = 0.8,
      text = "Proportion of haplotypes copied from all non-local donor ethnic groups\nby recipient East African Nilo-Saharan and Afroasiatic speakers")

#######################################
## scale bar

#tmp <- c(0,seq(0.1,max(copyhit4),length=10))
tmp <- seq(0,1,length=10)
par(mar=c(3,10,3,10))
image(0:length(tmp),
      1,
      as.matrix(tmp),
      col = heatcols,axes = F,
      xlab = "",ylab = "")
box()
xat <- pretty(0:length(tmp)) #c(0.5,2,4,6,8,1)
xlab <- pretty(tmp) # c(0,0.2,0.4,0.6,0.8,1)
axis(1,at = xat, labels = xlab, xpd = T)
mtext(3,cex = 0.8,
      text = "Prop of haplotypes copied from each donor group\n(normalised by number of inds in recipient group)")

#######################################################
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
mtext(2,at=ymax,text="c",adj=0, cex = 1.5, las = 2, line = 5)  

dev.off()

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


hlapops <- c(popkey$Ethnic_Group[popkey$RegionM=="Western_Africa_Niger-Congo"][1:10])
chrom <- "06"
plot_range <- c(30e6,35e6)
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
pdf("figures/Figure3HLA.pdf", height = 15, width = 14)
plot_matrix <- matrix(c(1:(2*(length(hlapops)+3))), 13,2,byrow = T)
layout(plot_matrix,
       heights = c(rep(1,length(hlapops)),1,1.5,2.5))

gene_region <- c(32.2e6,33e6)

for(pop in hlapops)
{
  print(pop)
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

  # chromlength <- as.numeric(as.character(map$recrate))
  # chromplot <- paintedchromreg[plot_region,]
  # chrompos <- as.numeric(as.character(map$position))[plot_region]
  # chromposI <- c(diff(chrompos),0)
  # chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
  # chromcols <- pcolshex[chromcolsbreaks]
  ## REORDER AROUND FOCAL SNP
  # focalsnp <- "rs9268831"
  # focalsnppos <- as.numeric(as.character(snps$pos[snps$rsid==focalsnp]))
  # focalsnppos <- which(chrompos==focalsnppos)
  # neword <- hclust(dist(t(chromplot[(focalsnppos-100):(focalsnppos+100),])))$order 
  # chromplot <- chromplot[,neword]
  # # image(chrompos,
  # #       1:ncol(chromplot),
  # #       (chromplot),
  # #       col=chromcols,breaks = c(0,chromcolsbreaks),
  # #       axes=F,xlab="",ylab="")
  # # mtext(3,line=0,text=" b",adj=0, cex = 2)
  # # 
  # tmp <- apply(chromplot,2,function(x){x==x[focalsnppos]})
  # for(i in 1:ncol(chromplot))
  # {
  #   print(i)
  #   switchs <-  which(diff(chromplot[focalsnppos:nrow(chromplot),i])!=0)[1]
  #   if(is.na(switchs)) switchs <- 1
  #   notcore <- (focalsnppos:nrow(chromplot))[switchs]:nrow(chromplot)
  #   tmp[notcore,i]   <- FALSE
  #   switchs <- which(diff(chromplot[1:focalsnppos,i])!=0)[1]
  #   if(is.na(switchs)) switchs <- focalsnppos
  #   notcore <- 1:switchs
  #   tmp[notcore,i]   <- FALSE
  # }
  # 
  # chromplot2 <- matrix(0,dim(chromplot)[1],dim(chromplot)[2])
  # chromplot2[tmp] <- chromplot[tmp]
  # neword <- hclust(dist(t(chromplot)))$order #[(focalsnppos-100):(focalsnppos+100),]
  # chromplot <- chromplot[,neword]
  # 
  # par(mar=c(0,6,1,1))
  # image(chrompos,
  #       1:ncol(chromplot),
  #       (chromplot2),
  #       col=c("white",chromcols),breaks = c(0,0.1,chromcolsbreaks),
  #       axes=F,xlab="",ylab="")
  # mtext(2,line=0,text=pop,adj=0.5, cex = 2, las = 2)
  # box()

  ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
  paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
  for(i in 1:nrow(paintedchromreg))
  {
    print(i)
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
  plot_region <- map$position>=plot_range[1]&map$position<=plot_range[2]
  chrompos <- map$position[plot_region]
  chromposI <- c(diff(chrompos),0)
  chromplot <- paintedchromregprop[plot_region,8:1]
  colnames(chromplot) <- rev(ancreg_list)
  chromplot <- chromplot/rowSums(chromplot)
  
  plot_reg <- "East_Africa_Nilo-Saharan"
  if(pop == "ANUAK") plot_reg <- "Western_Africa_Niger-Congo"
  y_max <- 0.35 #max(chromplot[,plot_reg])
  par(mar=c(0,6,0.5,1))
  plot(chrompos,chromplot[,1], type = "n",
       axes = F, xlab = "", ylab = "Ancestry Prop",
       ylim = c(0,y_max), xaxs = "i", yaxs = "i")
  
  abline(h=mean(chromplot[,plot_reg]), lty = 2)
  abline(v=gene_region, lty = 1, col = "grey", xpd =T)
  for(i in which(colnames(chromplot)%in%c(plot_reg,"Eurasia")))
  {
    points(chrompos,chromplot[,i],
         col = pcolshex[ancreg_list==colnames(chromplot)[i]],
         type = "S", lwd = 2)
  }
  xat <- pretty(chrompos)
  yat <- c(0,y_max)
  #if(pop == "ANUAK") axis(1, at = xat, labels = xat/1e6, xpd = T)  
  axis(2, at = yat, labels = yat, xpd = T, las = 2)
  box()
  legend("topright", legend = pop, bty = "n")
  par(mar=c(0,6,0.5,1))
  barplot(t(chromplot),
          width=chromposI,
          col=rev(pcolshex),xaxs="i",yaxs="i",
          space=0,axes=F,xaxt="n",border=NA,
          xlim=c(0,sum(chromposI)),horiz=F,
          xlab="")
  # xatlab <- pretty(chrompos)
  # xat <- seq(0,sum(chromposI), length.out = length(xatlab))
  # #xat[length(xat)] <- sum(chromposI)
  # axis(1,at=xat,labels=xatlab/1e6)
  # mtext(1,text=paste("position on chromosome",as.numeric(chrom), "(Mb)"),line=2,cex=1)
  # mtext(3,line=0,text=" d",adj=0, cex = 2)
  # mtext(4, text = "ancestry\nacross all\nhaplotypes",
  #       las = 2, adj = 0.5, line = 3.1, cex = 0.75)
  # 
  # 
}


## HAP 1257
## SHOW PAINTINGS ACROSS ALL POP OF INTEREST
psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
psamplesind <- (1:nrow(psamples))[psamples[,"ind"] == "ANUAK11"]
# 2 haps per sample!!
psampleshap <- hap2sampleindex(psamplesind,2)
psampleshap <- sort(c(psampleshap,psampleshap+1))
psamplesindsamp <- hap2sampleindex(psampleshap)
tmp <- c()
for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
paintedchrom <- datafile[paste0("/paintings/chrom",chrom,"/local")]
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
paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
for(i in 1:nrow(paintedchromreg))
{
  print(i)
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
# test_rev_region <- chrompos>=gene_region[1]&chrompos<=gene_region[2]
# tmp <- paintedchrom[test_rev_region,]
# 
####################################################################################
## SUBSAMPLE SNPS
#snpsample <- round(seq(1,nrow(chromplot), length = 1000))
####################################################################################
## PROPORTIONS    
chromplot <- paintedchromregprop[,8:1]
chromplot <- chromplot/rowSums(chromplot)
colnames(chromplot) <- rev(ancreg_list)
## ADJUST FOR RELEVANT PLOT REGION
plot_region <- chrompos>=plot_range[1]&chrompos<=plot_range[2]

par(mar=c(0,6,0.5,1))
plot(chrompos[plot_region],chromplot[plot_region,"Western_Africa_Niger-Congo"], 
     axes = F, xlab = "", ylab = "Ancestry Prop",
     ylim = c(0,1), xaxs = "i", yaxs = "i",
     type = "n", col = pcolshex[ancreg_list=="Western_Africa_Niger-Congo"])
abline(h=mean(chromplot[,plot_reg]), lty = 2)
abline(v=gene_region, lty = 1, col = "grey", xpd =T)
for(i in which(colnames(chromplot)%in%c("Western_Africa_Niger-Congo","Eurasia")))
{
  points(chrompos[plot_region],chromplot[plot_region,i],
         col = pcolshex[ancreg_list==colnames(chromplot)[i]],
         type = "S", lwd = 2)
}
xat <- pretty(chrompos)
yat <- c(0,1)
#if(pop == "ANUAK") axis(1, at = xat, labels = xat/1e6, xpd = T)  
axis(2, at = yat, labels = yat, xpd = T, las = 2)
box()


par(mar=c(0,6,0.5,1))
barplot(t(chromplot[plot_region,]),
        width=chromposI[plot_region],
        col=rev(pcolshex),xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        horiz=F,#xlim=c(0,sum(chromposI)),
        xlab="", ylab = paste0("ANUAK 11\nancestry"),
        xaxs = "i", yaxs = "i")

xat <- pretty(c(0,sum(chromposI[plot_region])))
xlab <- pretty(chrompos[plot_region])/1e6
axis(1, at = xat[1:length(xlab)], labels = xlab)
#######################################################
## GENES and RECRATES
#######################################################
### Code to load genetic maps from Gav
library(hapdb)
library(RSQLite)
filename = "/mnt/kwiat/data/1/galton/malariagen/human/reference/genetic_maps/maps_b37/maps_b37.sqlite"
genetic.map.db = dbConnect( dbDriver( "SQLite" ), filename )
gm = dbGetQuery(
  genetic.map.db,
  sprintf(
    "SELECT * FROM GeneticMap WHERE chromosome == '%s' AND position >= %d AND position <= %d ORDER BY position",
    chrom, plot_range[1], plot_range[2]
  )
)

gm$COMBINED_LD_rate = NA
gm$COMBINED_LD_rate[-1] = ( gm$COMBINED_LD[-1] - gm$COMBINED_LD[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )
gm$YRI_LD_rate = NA
gm$YRI_LD_rate[-1] = ( gm$YRI_LD[-1] - gm$YRI_LD[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )
gm$African_Enriched_rate = NA
gm$African_Enriched_rate[-1] = ( gm$African_Enriched[-1] - gm$African_Enriched[-nrow(gm)] ) / ( gm$position[-1] - gm$position[-nrow(gm)] )

for(plot_run in 1:2)
{
  par(mar=c(0,6,1,1))
  plot(gm$position[gm$chromosome==chrom],gm$COMBINED_LD_rate[gm$chromosome==chrom]*1e6,
       xlim=range(plot_range), type = "S", ylim = c(0,100),
       axes = F, ylab = "Rec Rate\n(cM/Mb)", xlab = "", xaxs = "i")
  axis(2, las = 2)
  abline(v=gene_region, col = "grey", xpd = T)
  
  points(gm$position[gm$chromosome==chrom],gm$YRI_LD_rate[gm$chromosome==chrom]*1e6,
         type = "S", col = "blue", lty = 2)
  points(gm$position[gm$chromosome==chrom],gm$African_Enriched_rate[gm$chromosome==chrom]*1e6,
         type = "S", col = "red", lty = 3)
  legend("topleft",legend=c("HAPMAP","YRI", "Africa enriched"),
         col = c("black","blue","red"), lty = c(1,2,3), bty = "n", horiz = T)

}
#mtext(2,at=100,text=panel_letter,adj=0, cex = 1.5, las = 2, line = 5)          
#######################################################
### PLOT GENES
source("~/repos/glycophorins/external_software/plot_genes.R")
genes = load.genes( "/mnt/kwiat/data/1/galton/malariagen/human/reference/genome-mysql.cse.ucsc.edu/2015-08-18/UCSC_hg19_2015-08-18_refGene.tsv" )
for(plot_run in 1:2)
{
  par(mar=c(4,6,1,1))
  plot.genes( chromosome = chrom, region = range(pretty(gene_region)), genes, xaxt = "n" , plot.ylab = "Genes")
  x_at <- pretty(gene_region)
  axis(1,at=x_at,labels=x_at/1e6, xpd = T)
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
dev.off() 
  
####################################################################  
####################################################################
####################################################################  
  ## WHO'S BEING COPIED??
  tmp <- paintedchrom[chrompos>=plot_range[1]&chrompos<=plot_range[2],]
  happrops <- matrix(0,nr=nrow(tmp), ncol = length(happops))
  colnames(happrops) <- 1:ncol(happrops)
  for(i in 1:nrow(tmp))
  {
    tmp2 <- table(tmp[i,])
    happrops[i,match(names(tmp2),colnames(happrops))] <- tmp2
  }
  
  chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
  chromcols <- pcolshex[chromcolsbreaks]
  chromplot <- paintedchromreg

  
### HAP 1257
## SHOW PAINTINGS ACROSS ALL POP OF INTEREST
psamplesind <- (1:nrow(psamples))[psamples[,"region"] == pop]
psamplesind <- (1:nrow(psamples))[psamples[,"ind"] == "ANUAK11"]
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
paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
for(i in 1:nrow(paintedchromreg))
{
  print(i)
  tmp <- table(paintedchromreg[i,])
  paintedchromregprop[i,as.numeric(names(tmp))] <- tmp
}
  
chromlength <- as.numeric(as.character(map$recrate))
chrompos <- as.numeric(as.character(map$position))
chromposI <- c(diff(chrompos),0)
chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
chromcols <- pcolshex[chromcolsbreaks]
chromplot <- paintedchromreg

  
## WHO'S BEING COPIED??
tmp <- paintedchrompop[chrompos>=plot_range[1]&chrompos<=plot_range[2],]
  popprops <- matrix(0,nr=nrow(tmp), ncol = length(popkey$Ethnic_Group))
  colnames(popprops)<- popkey$Ethnic_Group
  for(i in 1:nrow(tmp))
  {
    tmp2 <- table(tmp[i,])
    popprops[i,match(names(tmp2),colnames(popprops))] <- tmp2
  }
  
  popprops <- popprops[,colnames(popprops)%in%popkey$Ethnic_Group[popkey$RegionM == copy_reg]]
  
  hlacopy <- rbind(hlacopy,cbind(pop,popprops))
  
  
 
####################################################################################
# test_rev_region <- chrompos>=gene_region[1]&chrompos<=gene_region[2]
# tmp <- paintedchrom[test_rev_region,]
# 
####################################################################################
## SUBSAMPLE SNPS
#snpsample <- round(seq(1,nrow(chromplot), length = 1000))
####################################################################################
## PROPORTIONS    
chromplot <- paintedchromregprop[,8:1]
chromplot <- paintedchromreg
## ADJUST FOR RELEVANT PLOT REGION
plot_region <- chrompos>=plot_range[1]&chrompos<=plot_range[2]

par(mar=c(0.5,6,0,1))
barplot(t(chromplot[plot_region,]),
        width=chromposI[plot_region],
        col=rev(pcolshex),xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        horiz=F,#xlim=c(0,sum(chromposI)),
        xlab="", ylab = paste0(pop,"\nancestry"))
# abline(v=(seq(0,sum(chromposI[plot_region]),
#               length=sum(plot_region)))[chrompos[plot_region]==gene_region[1]],
#           col = "grey10", lty = 1, xpd = T)
# abline(v=(seq(0,sum(chromposI[plot_region]),
#               length=sum(plot_region)))[chrompos[plot_region]==gene_region[2]],
#        col = "grey10", lty = 1, xpd = T)
# #######################################################

## add a plot on top to add lines
# par(new = T)
# plot(0,0,axes = F, ylim=c(0,1), type = "n",
#      xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
# abline(v=gene_region, col = "grey10", lty = 1, xpd = T)


y_max <- max(rowSums(popprops))
x_pos <- chrompos[chrompos>=31.5e6&chrompos<=33.0e6]


par(mar=c(4,4,2,2))
plot(x_pos,rowSums(popprops),
     ylim=c(0,y_max), type = "S", axes = F,
     xlab = "pos", ylab = "prob", lwd = 2)
axis(1)
axis(2, las = 2)
for(i in 1:ncol(popprops))
{
  points(x_pos,popprops[,i], col = c("red","blue","green","hotpink")[i], type = "S", lwd = 2)
}





#

clengths <- c()
for(pop in hlapops)
{
  pop <- "ANUAK"
  print(pop)
  in_file <- paste0("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/paintingstats/",pop,"Chrom06paintingstatsByReg.gz")
  tmp <- read.table(in_file,header = T, as.is = T, stringsAsFactors = F)
  clengths <- rbind(clengths,tmp)
}

clengths$length <- clengths$endpos-clengths$startpos
xmax <- 2e6
xbreaks <- seq(0,max(clengths$length), by = 2.5e4)
hist(clengths$length[clengths$poppainting=="Western_Africa_Niger-Congo"],
     col = "red", xlim = c(0,xmax), breaks = xbreaks,
     axes = F, xlab = "lengths copied from ANUAK (Mb)",
     xaxs = "i", yaxs = "i")
x_at <- pretty(c(0,xmax))
axis(1,at = x_at, labels = x_at/1e6)
axis(2, las = 2)
## ADD HLA
hist(clengths$length[clengths$startpos>=32e6&clengths$endpos<=33e6&clengths$poppainting=="Western_Africa_Niger-Congo" ],
     breaks = xbreaks, col = "blue", xlim = c(0,2e6), add = T)


###
### SCRIPT TO PLOT LAMBDAS and PVALS EQUIVALENT TO MVN ANCESTRY PLOTS ###
###
res_dir <- "/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/"
#res_dir <- "/well/malariagen/malariagen/human/george/copy_selection/"
source("~/R/Copy/Rprojects/AfricaPOPGEN/functions/plottingFunctions.R")
#source("plottingFunctions.R")
source("/mnt/kwiat/well/human/george/copy_selection/R/hitplots.R")
#source("hitplots.R")
snps <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330Kphased.legend", header=F,as.is=T)
colnames(snps) <- c("chrom","rsid","pos","a0","a1")

info3 <- read.table("~/R/Copy/Rprojects/AfricaPOPGEN/admixture/MalariaGenAdmixturePopulationKey3.txt",header=T,comment.char="")
#info3 <- read.table("MalariaGenAdmixturePopulationKey3.txt",header=T,comment.char="")
###################
## COLOUR VECTOR ##
pcols <- c("darkgreen","darkolivegreen1","darkorchid","navyblue","red","hotpink","cadetblue1")
pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00")
pop_dir <- "/mnt/kwiat/well/human/george/chromopainter2/analysislists/"
#pop_dir <- "/well/malariagen/malariagen/human/george/chromopainter2/analysislists/"
popkey <- read.table(paste0(pop_dir,"populationOverviewMalder.txt"),header=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"
regions <- as.character(unique(popkey$RegionM))
minp <- 20
pop <- "FULAI"

in_root <- paste0(res_dir,pop,"nolocalChrom")
in_suff <- ".ancestryselectionVII.gz"
pl <-c()
for(chrom in 1:22)
{
  if(chrom < 10) chrom <- paste0("0",chrom)
  in_file <- paste0(in_root,chrom,in_suff)
  if(chrom == 1)
  {
      pl <- read.csv(in_file,header=T, as.is=T)
  } else
  {
      pl <- rbind(pl,read.csv(in_file,header=T, as.is=T))
  }
}


########################################################################
source("~/R/Copy/Rprojects/AfricaPOPGEN/functions/hitplots.R")
genes <- load.genes()
############################################################
## TO DO:

## GENERATE PVAL USING LRT FOR EACH ANCESTRY AT EACH SNP
## USE CODE FROM MVN PLOT TO MAKE SAME PLOT, BUT
##      PLOT 6 X P VALUES -- WITH DIFFERENT COLOURED LINES
##      INSTEAD OF ANCESTRY DIFFERENCES, PLOT LAMBDAS AT EACH SNP?


## PVALS
plotregs <- gsub(".GB.beta","",colnames(pl)[grep(".GB.beta",colnames(pl))])
n_plotregs <- length(plotregs)
pval <- c()
lams <- c()
## get lambdas
ls <- pl[,grep(".GB.beta",colnames(pl))]
ps <- pl[,grep(".GB.P",colnames(pl))]
colnames(ls) <- colnames(ps) <- paste(pop,gsub(".GB.beta","",colnames(ls)),sep=".")
if(is.null(lams))
{
    lams <- ls
    pval <- ps
} else 
{
    lams <- cbind(lams,ls)
    pval <- cbind(pval,ps)
}
## some plotting parameters ##
regname <- "LCT"
grepgene <- c(regname,"MCM6")
chrom <- 2
lower_pos <- 132e6
higher_pos <- 142e6
plot_pops <- "FULAI"
#pops <- c("CEU","GBR","FIN","TSI","IBS","CHB","GIH")
pops <- "FULAI"
ymax <- 20
rightmar <- 0

# ##Define region to plot
regname <- "DARC"
grepgene <- c(regname,"KIRREL")
chrom <- 1
lower_pos <- 153e6
higher_pos <- 163e6
plot_pops <- "FULAI"
#pops <- c("CEU","GBR","FIN","TSI","IBS","CHB","GIH")
pops <- "FULAI"
ymax <- 20
# 


#outpdf  <-paste0("R/figures/",pop,"logRatioTestGenes",regname,".pdf")
#pdf(outpdf,width=12,height=12)

outpng <- paste0("R/figures/",pop,"logRatioTestGenes",regname,"2.png")
png(outpng,width=3000,height=3000,res=300)
npops <- 1# length(pops)
layout(matrix(c(3,6,4,5,1,7,2,8),byrow=T,4,2),
       heights=c(1.5,5,5,0.5),
       widths=c(12,0))

## plot region
test<-snps$chrom==chrom&snps$pos>lower_pos&snps$pos<higher_pos

## plot pvals
par(mar=c(2,4,0,rightmar))
for(i in 1:ncol(pval))
{
    reg <- gsub(paste0(pop,"."),"",colnames(pval)[i])
    reg <- gsub("\\.","-",reg)
    regcol <- pcolshex[regions==reg]
    if(i ==1)
    {
        plot(snps$pos[test],pval[test,i],
             ylab="",xlab="",axes=F,
             type="S",lwd=2,ylim=c(0,ymax),col=regcol)
    } else
    {
        points(snps$pos[test],pval[test,i],type="S",lwd=3,col=regcol)
    }
}
axis(2,las=2,cex.axis=3,line=-3)
mtext(2,text="-log10 P",line=1,cex=2)
axis(1,at=x$pos[test],labels=NA,tck=0.05,line=1,lwd=0,lwd.ticks=1)
abline(h=7,lty=2,xpd=F)

##Plot genes
#Genes
chrom2 <- chrom
if(chrom2<10) chrom2 <- paste("0",chrom,sep="")
genes1 <- genes[ which( genes$chromosome == chrom2 & genes$txEnd >= lower_pos & genes$txStart <= higher_pos ),, drop = FALSE ]
genes2grep <- c()
for(i in 1:length(grepgene)) genes2grep <- c(genes2grep,grep(grepgene[i],genes1$name2))
genes1 <- genes1[ genes2grep,]

par(mar=c(0,4,0,rightmar))
#plot.genes(c(lower_pos,higher_pos), genes1, 12,set.gene.cex=T,gene.max=3)
plot.genes(c(lower_pos,higher_pos), genes1, 2,set.gene.cex=T,gene.max=2,gene.cex=2)
mtext(2,text="select genes",adj=0,line=1,cex=2,las=1)

## top legedn
par(mar=c(0,0,0,0))
plot(0,0,type="n",axes=F,xlab="",ylab="")
par( lend = 1 ) 
plot_reg <- sort(match(gsub("\\.","\\-",gsub(paste0(pop,"."),"",colnames(lams))),regions))
legend("bottom",
       legend=gsub("\\_"," ",regions[plot_reg]),
       ncol=2,cex=2.5,pt.cex=4,bty="n",
       col=pcolshex[plot_reg],
       lwd=4,lty=1)

## PLOT LAMBDAS
par(mar=c(4,4,0,rightmar))
ymins <- max(abs(lams))
for(i in 1:ncol(lams))
{
    reg <- gsub(paste0(pop,"."),"",colnames(lams)[i])
    reg <- gsub("\\.","-",reg)
    regcol <- pcolshex[regions==reg]
    if(i ==1)
    {
        plot(snps$pos[test],lams[test,i],
             ylab="",xlab="",axes=F,
             type="S",lwd=2,ylim=c(-ymins,ymins),col=regcol)
    } else
    {
        points(snps$pos[test],lams[test,i],type="S",lwd=3,col=regcol)
    }
}
axis(2,las=2,cex.axis=3,line=-3)
mtext(2,text=expression(Delta ~ ancestry),line=1,cex=2)
abline(h=0,lty=2,xpd=F)
axis_seq <- seq(lower_pos,higher_pos,length=6)
axis(1,las=2,at=axis_seq,labels=rep("",length(axis_seq)),las=1,cex.axis=3,lwd=1,line=0)
axis(1,las=2,at=axis_seq,labels=axis_seq/1e6,las=1,cex.axis=3,lwd=0,line=1)
mtext(1,text=paste("position on chromosome",chrom,"(Mb)"),cex=2,line=4)

dev.off()





###############################################
## CODE TO LOOK AT ANCESTRY SELECTION OUTPUT ##
###############################################

library(dplyr) ## for binding dataframes with difference numbers of columns; bind_rows
popkey <- read.table("data/PopulationKey.txt",
                     header = T, stringsAsFactors = F, as.is= T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
## going to work with chrom 22 from the FULAI for the time being
indir <- "/mnt/kwiat/data/1/galton/users/george/copy_selection/output/"
pop <- "FULAI"
chrom <- 10
in_file <- paste0(indir,pop,"nolocalChrom",chrom,".ancestryselectionIII.gz")

lrts <- read.csv(in_file, header = T)

popstab <- read.table("data/PopulationKey.txt", header = T)
pops <- toupper(as.character(levels(popstab$Ethnic_Group)))
regions <- as.character(levels(popstab$AncestryRegion))
pcols_table <- read.table("~/repos/ANCSELECT/data/RegionColours.txt",
                          comment.char = "", header = T)


# tmp <- read.csv("/mnt/kwiat/well/human/george/copy_selection/Copy_Selection_deviant_FULAI.csv", header = T)
# tmp <- tail(tmp,3392)
# mvn <- -log10(tmp$pval)
mvn <- -log10(lrts$MVNp)

y_lims <- c(0,15)
y_at <- pretty(y_lims)
x_at <- pretty(lrts$pos)
x_at[which.min(x_at)] <- min(lrts$pos)
x_at[which.max(x_at)] <- max(lrts$pos)
x_lims <- range(x_at)



#########################################################


n_plots <- length(regions)+2

png("figures/ModelComparisonFULAIChrom22.png",
    width = 2400,height = 200*(n_plots), res = 300)

plot_matrix <- matrix(cbind(c(n_plots,1:(n_plots-1)),
                            c((n_plots+1):(n_plots+10))),
                      nc = 2)
layout(plot_matrix,widths = c(10,2)) 

plot_x <- mvn
plotcol <- "black"

par(mar = c(1,4,0,1))
plot(lrts$pos,plot_x,pch = 20, col = plotcol, type = "S",
     ylim = y_lims, ylab = "-log10 P", xlab = "",
     axes = F, xlim = x_lims, xaxs = "i", yaxs = "i")
axis(2, las = 2, at = y_at, labels = y_at)
legend("topright", bty = "n", legend = "combined MVN", text.col = plotcol)

for(reg in regions)
{
  reg_col <- paste0(gsub("\\-","\\.",reg),".P")
  if(sum(colnames(lrts)==reg_col)>0)
  {
    plot_x <- lrts[,reg_col]
    plotcol <- as.character(pcols_table$Colour[pcols_table$AncestryRegion==reg])
    plot(lrts$pos,plot_x,pch = 20, col = plotcol, type = "S",
         ylim = y_lims, ylab = "-log10 P", xlab = "",
         axes = F, xlim = x_lims, xaxs = "i", yaxs = "i")
    axis(2, las = 2, at = y_at, labels = y_at)
    legend("topright", bty = "n", legend = gsub("\\_"," ",reg), text.col = plotcol)
  }
}

par(mar = c(5,4,0,1))
plot(lrts$pos,plot_x,pch = 20, col = plotcol, type = "n",
     ylim = y_lims, ylab = "", xlab = "",
     axes = F, xlim = x_lims, xaxs = "i", yaxs = "i")
axis(1, at = x_at, labels = round(x_at/1e6), xpd = T)
mtext(1,line=2,text = paste0("position on chromosome ",chrom," (Mb)"))


#########################################################
## 01 PLOT A PAINTED CHROMOSOME
library("h5")
## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
  ## finds the first sample index for a haplotype
  sample <- (hap*nsamps)-(nsamps-1)
  return(sample)
}

popkey_file <- "~/repos/popgen/data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)

pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )


#########################################################
## 01 PLOT A PAINTED CHROMOSOME
plot_range <- x_lims
## GET MAP AND POSITION INFO
map <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/map")]))
colnames(map) <- c("position","recrate")
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- readDataSet(datafile[paste0("/paintings/samples/individuals")])
colnames(psamples) <- c("ind","region","X")

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
## ADJUST FOR RELEVANT PLOT REGION
plot_region <- chrompos>=plot_range[1]&chrompos<=plot_range[2]

par(mar=c(0.5,4,0.5,1))
bp <- barplot(t(chromplot[plot_region,]),
        width=chromposI[plot_region],
        col=rev(pcolshex),xaxs="i",yaxs="i",
        space=0,axes=F,xaxt="n",border=NA,
        horiz=F,# xlim=c(0,sum(chromposI)),
        xlab="", ylab = paste0(pop,"\nancestry"))

###################################################################################
## COMPARE MVN MARGINAL PVALUES
plot(0,0, type = "n", xlab = "", ylab = "", axes = F)
par(mar = c(0,0,0,0))
plot(0,0, type = "n", xlab = "", ylab = "", axes = F)
legend("top", legend=c("LRT v MVN pvalue", "LRT v marginal MVN pvalue"),
       col = "black", bty= "n",pch = c(22,21), cex = 0.75)

for(reg in regions)
{
  reg_col <- paste0(gsub("\\-","\\.",reg),".P")
  if(sum(colnames(lrts)==reg_col)>0)
  {
    
    par(mar = c(1,1,0,1))
    plot(0,0, type = "n", xlim=y_lims,ylim= y_lims,axes = F,
         ylab = "", xlab = "",
         xaxs = "i", yaxs = "i")
    abline(a=0,b=1, col = "red")
    axis(4, at = y_at, las = 2, tck = +0.025, labels = NA)
    axis(1, las = 1, tck = +0.025, labels = NA)
    box()
    if(reg == "South_Africa_Niger-Congo") axis(1)
    
    reg_col2 <- paste0(gsub("\\-","\\.",reg),".MVNp")
    reg_col3 <- "MVNp"
    for(plot_y in c(reg_col2,reg_col3))
    {
      plot_pch <- 21
      if(plot_y == reg_col3) plot_pch <- 22
      plot_x <- lrts[,reg_col]
      plot_y <- -log10(lrts[,plot_y])
      plotcol <- as.character(pcols_table$Colour[pcols_table$AncestryRegion==reg])
      points(plot_y,plot_x,pch = plot_pch, bg = plotcol)
    }
  }
}


dev.off()

png("figures/ModelComparisonFULAIChrom22RyanGeorgeLRT.png",
    width = 750, height = 1500, res = 150)
layout(matrix(c(1:(n_plots-2)),nc = 2))
for(reg in regions)
{
  reg_col <- paste0(gsub("\\-","\\.",reg),".P")
  if(sum(colnames(lrts)==reg_col)>0)
  {
    par(mar = c(4,4,1,2))
    plot(0,0, type = "n", xlim=c(0,8),ylim=c(0,8),axes = F,
         ylab = "Ryan's LRT", xlab = "LRT",
         xaxs = "i", yaxs = "i")
    abline(a=0,b=1, col = "red", lwd = 2)
    box()
    axis(1)
    axis(2, las = 2)
    
    reg_col2 <- paste0(gsub("\\-","\\.",reg),".PII")
    plot_pch <- 21
    plot_x <- lrts[,reg_col]
    plot_y <- lrts[,reg_col2]/10
    plotcol <- as.character(pcols_table$Colour[pcols_table$AncestryRegion==reg])
    points(plot_x,plot_y,pch = plot_pch, bg = plotcol, xpd = T)
    mtext(4, text = gsub("\\_"," ",reg), col = plotcol, cex = 0.5)
  }
}

dev.off()


png("figures/ModelComparisonFULAIChrom22LRTqqPLOTS.png",
    width = 750, height = 1500, res = 150)
layout(matrix(c(1:(n_plots-2)),nc = 2))
for(reg in regions)
{
  reg_col <- paste0(gsub("\\-","\\.",reg),".P")
  if(sum(colnames(lrts)==reg_col)>0)
  {
    par(mar = c(4,4,1,2))
    plot_pch <- 21
    plot_x <- lrts[,reg_col]
    plotcol <- as.character(pcols_table$Colour[pcols_table$AncestryRegion==reg])
    ax_lims <- c(0,max(plot_x))
    plot(sort(-log10(ppoints(n_snps))),
         sort(plot_x),
         pch = plot_pch, bg = plotcol, xpd = T,
         ylab = "LRT pvalues", xlim = ax_lims,ylim =ax_lims,
         xlab = "Theoretical pvalues", yaxt = "n")
    axis(2, las = 2)
    abline(a=0,b=1,col = "red", lwd = 2)
    
    mtext(4, text = gsub("\\_"," ",reg), col = plotcol, cex = 0.5)
  }
}

dev.off()





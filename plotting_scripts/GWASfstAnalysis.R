#######################################################
### SCRIPT TO PLOT FST ACROSS GWAS ETHNIC GROUPS    ###
#######################################################
library("dplyr")

#############################################
## GET SNPS
snps <- read.table("/mnt/kwiat/well/human/george/GWAS_painting/data/Superset2_GWAS-2.5M.snps",
                   header = F, as.is = T)
colnames(snps) <- c("chrom","snp","rsid","pos","a0","a1")
#############################################
## GET FSTS
indir <- "/mnt/kwiat/well/human/george/GWAS_painting/fst/"
pops <- c("Malawi","Kenya","Gambia")
ethnics <- c("MALAWI","CHONYI","GIRIAMA","KAUMA","FULA","WOLLOF","MANDINKA","JOLA")
pcols <- c("brown","darkorange","darkgoldenrod1", "coral","darkgreen","cornflowerblue","darkblue","cadetblue")
chroms <- 1:22

fsts <- matrix(0,nc=0,nr=nrow(snps))
for(pop1 in ethnics)
{
  for(pop2 in ethnics)
  {
    if(pop1 != pop2 & !paste(paste(pop2,pop1,sep="-"),"num",sep = ".") %in% colnames(fsts))
    {
      fs <- c()
      for(chrom in chroms)
      {
        if(chrom < 10) chrom <- paste("0",chrom, sep = "")
        infile <- paste(indir,pop1,"_",pop2,"_Chrom",chrom,"hudsonfst.txt", sep = "")
        if(file.exists(infile))
        {
          tmp <- read.table(infile, header = T, as.is = T, stringsAsFactors = F)  
          fs <- rbind(fs,tmp)
        } else
        {
          infile <- paste(indir,pop2,"_",pop1,"_Chrom",chrom,"hudsonfst.txt", sep = "")
          if(file.exists(infile))
          {
            tmp <- read.table(infile, header = T, as.is = T, stringsAsFactors = F)
            ## SWITCH SO THAT THE COLUMNS ARE IN THE CORRECT ORDER
            pop1mat <- tmp[,c("a0.pop2","a1.pop2")]
            pop2mat <- tmp[,c("a0.pop1","a1.pop1")]
            tmp[,c("a0.pop1","a1.pop1")] <- pop1mat
            tmp[,c("a0.pop2","a1.pop2")] <- pop2mat
            fs <- rbind(fs,tmp)
          }
        }
      }
      print(paste(paste(pop1,pop2,sep="_"),nrow(fs)))
      ## NOW CHANGE COLUMN HEADERS AND ADD TO MASTER FST MATRIX
      colnames(fs) <- gsub("pop1",pop1,colnames(fs))
      colnames(fs) <- gsub("pop2",pop2,colnames(fs))
      for(hd in c("num","den","hudson.fst"))
      {
        colnames(fs) <- gsub(hd,paste(paste(pop1,pop2,sep="-"),hd,sep = "."),colnames(fs))
      }
      if(nrow(fs) == nrow(snps))
      {
        fsts <- cbind(fsts,fs[,!colnames(fs)%in%colnames(fsts)])
      }
    }
  }
}

fsts.back <- fsts
## CONVERT TO dplyr DATA TABLE
fsts <- tbl_df(fsts)

## LOOK AT THE DATA
glimpse(fsts)

#############################################
## QUICK SEQUE INTO MAKING FILES FOR TREESELECT
## EACH FILE HAS FOR A GIVEN POPULATION:
##  rsid, chrom, position, N, A1, A2, A1-freq
# tsdir <- '/mnt/kwiat/well/human/george/GWAS_painting/treeselect/data/'
# for(pop in ethnics)
# {
#   N <- fsts[,paste0("a0.",pop)] + fsts[,paste0("a1.",pop)]
#   a0freq <- fsts[,paste0("a0.",pop)] / N
#   tmp <- cbind(snps[,c("rsid","chrom","pos")],N,snps[,c("a0","a1")], a0freq)
#   colnames(tmp) <- c("SNP-id","CHR","POS","N","A1","A2","A1-frequency")
#   write.table(tmp,file=paste0(tsdir,pop,"treeselect.txt"),
#               sep = "\t",col.names = T, row.names = F, quote = F)
# }
### AND USE GAVIN'S CODE TO GET THE DATA FOR EUROPEAN 1KGP POPS
# library(RSQLite)
# allelesfile <- "/mnt/kwiat/data/1/galton/malariagen/human/data/callsets/sequenced/MalariaGEN_1000GP_combined_reference_panel_v1b/stats/allele_counts.sqlite"
# db <-  dbConnect( dbDriver( "SQLite" ), allelesfile )
# var1 <- snps$rsid
# 
#counts <- dbGetQuery( db, paste("SELECT * FROM ByPOPView LIMIT 10") ) ;

#############################################
## PRELIMINARY PLOTS
n_plots <- length(ethnics)
y_max <- 0.3 #round(max(fsts),2)
################################
## SORT OUT AN X-AXIS SCALE ###
## THIS IS A HACK FROM STEPHEN TURNERS MANHATTAN PACKAGE ##
pos2 <- rep(NA,nrow(snps))
labs <- unique(snps$chrom)
nchr <- length(chroms)
lastbase <- 0
ticks <- NULL
for (i in labs)
{
  if (i == 1)
  {
    pos2[snps$chrom == i] <- snps$pos[snps$chrom == i]
  } else
  {
    lastbase <- lastbase + tail(snps$pos[snps$chrom==(i-1)],1)
    pos2[snps$chrom == i] <- snps$pos[snps$chrom == i] + lastbase
  }
  ticks = c(ticks, pos2[snps$chrom==i][floor(length(pos2[snps$chrom==i])/2) + 1])
}

xmax <- ceiling(max(pos2) * 1.03)
xmin <- floor(max(pos2) * -0.03)
################################


png("figures/GWASpaintingFst.png",height = 2000, width = 2000, res = 150)
layout(matrix(1:n_plots))
for(j in ethnics)
{
  print(paste0("plotting ", j))
  par(mar=c(1,6,0.5,2))
  plot(pos2,rep(0,nrow(snps)),
       ylab="",xlab="",
       xaxt="n",type="n",yaxt="n",
       ylim=c(0,y_max),xlim=c(xmin,xmax))
  mtext(2,text=expression(italic(F[ST])), las = 2, line= 2.5)
  mtext(4,las=0,text=j,line=-2, col = pcols[ethnics == j])
  for(i in unique(labs)) abline(v=tail(pos2[snps$chrom<=i],1),lty=2)
  axis(2,las=2)  
  mtext(3,text=chroms,at=ticks,  line = -2, col = "grey")
  
  ## ORDER COLUMNS THAT WE WANT TO PLOT IN REVERSE
  select(fsts, contains(j)) %>%
    select(contains("hudson.fst")) %>%
    apply(2,mean, na.rm=T) %>%
    sort(decreasing=T) -> plot_order

  for(k in names(plot_order))
  {
    pop2 <- gsub(".hudson.fst","",gsub(paste0("-",j),"",gsub(paste0(j,"-"),"",k)))
    pop2col <- pcols[ethnics == pop2]  
    ploty <- select(fsts,matches(k))
    ## plot only top 5 percent
    plot_cex <- rep(1,nrow(ploty))
    #plot_cex[t(ploty)<=quantile(t(ploty),0.999, na.rm = T)] <- 0
    points(pos2,y=t(ploty),col=pop2col, pch = 20, cex = plot_cex) #type="S" #as.vector(t(ypnts))
    #points(pos2,y=t(ploty),col=pop2col, pch = 20, type="S") #as.vector(t(ypnts))
  }
}

dev.off()

##################################
## CODE BELOW ADAPTED FROM:
## http://personal.tcu.edu/kylewalker/interactive-flow-visualization-in-r.html?utm_content=bufferba958&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
## FIRST MAKE NEW DATAFRAME WITH PAIRWISE FST COMPS
plot_parset <- F
if(plot_parset == T)
{
  poppairfsts <- c()
  for(pop1 in ethnics)
  {
    for(pop2 in ethnics[ethnics!=pop1])
    {
      popcomb <- paste0(pop1,"-",pop2)
      colmnname <- paste0(popcomb,".hudson.fst")
      if(!colmnname%in%colnames(fsts)) popcomb <- paste0(pop2,"-",pop1)
      popfst <- apply(select(fsts,contains(popcomb)),2,sum, na.rm=T)
      popfst <- popfst[1]/popfst[2]
      poppairfsts <- rbind(poppairfsts,c(pop1,pop2,popfst))
    }
  }
  colnames(poppairfsts) <- c("pop1","pop2","fst")
  poppairfsts$fst <- as.numeric(poppairfsts$fst)
  poppairfsts <- tbl_df(poppairfsts)
  
  poppairfsts$pop1 <- factor(poppairfsts$pop1,levels=ethnics)
  poppairfsts$pop2 <- factor(poppairfsts$pop2,levels=rev(ethnics))
  poppairfsts <- poppairfsts[order(poppairfsts$pop1,poppairfsts$pop2),]
  
  library(parsetR) # devtools::install_github("timelyportfolio/parsetR")
  p <- parset(poppairfsts, dimensions = c('pop1', 'pop2'), 
         value = htmlwidgets::JS("function(d){return d.fst}"), 
         tension = 0.5, spacing = 50)
  p
}


#################################
## NOW COMPARE WITH ABACUS BAYESFACTORS

country <- "Gambia"
pop1 <- "FULA"
indir <- paste0("/mnt/kwiat/well/human/george/GWAS_painting/abacus/",country,"/",pop1,"/")
bayesf <- c()
for(chrom in 1:7) #chroms)
{
  if(chrom < 10) chrom <- paste0("0",chrom)
  tmp <- read.table(paste0(indir,chrom,"/",pop1,"GWASphasedChrom",chrom,".bayesfactor"),
                    header = F, as.is = T, stringsAsFactors = F)
  bayesf <- rbind(bayesf,tmp)
}
colnames(bayesf) <- c("snp","pos","bf","odds")

y_max <- 20

png(paste0("figures/GWASpaintingABACUS",pop1,".png"),
           height = 2000/8, width = 2000, res = 150)
  par(mar=c(1,6,0.5,2))
  plot(pos2,rep(0,nrow(snps)),
       ylab="",xlab="",
       xaxt="n",type="n",yaxt="n",
       ylim=c(0,y_max),xlim=c(xmin,xmax))
  mtext(2,text=expression(italic(-log[10]~BF)), las = 2, line= 2.5)
  mtext(4,las=0,text=pop1,line=-2, col = pcols[ethnics == pop1])
  for(i in unique(labs)) abline(v=tail(pos2[snps$chrom<=i],1),lty=2)
  axis(2,las=2)  
  mtext(3,text=chroms,at=ticks,  line = -2, col = "grey")
  points(pos2[1:length(bayesf$bf)],y=bayesf$bf,type="S",col=pcols[ethnics == pop1]) #as.vector(t(ypnts))
dev.off()


###########################################
## COMPUTE WINDOWED STATISTIC
## FOR HUDSON'S
##  - WINDOW THE GENOME BY DISTANCE
##  - THEN sum(num)/sum(den) OVER WINDOW
##  - COMPUTE JACKNIFE FOR SE?
## FOR ABACUS
##  - ASK GAVIN WHAT'S THE BEST WAY TO
##    AVERAGE BAYES FACTORS ACROSS A WINDOW







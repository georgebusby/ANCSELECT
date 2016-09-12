#######################################################
### SCRIPT TO WORK WITH VITERBI LENGTHS AT SNPS     ###
#######################################################

library(ggplot2)
library(dplyr)


#ethnics <- c("MALAWI","CHONYI","GIRIAMA","KAUMA","FULA","WOLLOF","MANDINKA","JOLA")
#pcols <- c("brown","darkorange","darkgoldenrod1", "coral","darkgreen","cornflowerblue","darkblue","cadetblue")

ethnics <- c("GIRIAMA","KAUMA","FULA","WOLLOF","MANDINKA","JOLA")
countries <- c("Kenya","Kenya","Gambia","Gambia","Gambia","Gambia")
pcols <- c("darkorange","coral","darkgreen","cornflowerblue","darkblue","cadetblue")

indir <- '/mnt/kwiat/well/human/george/GWAS_painting/length_analysis/allele_lengths/'
ancdir <- '/mnt/kwiat/well/human/george/GWAS_painting/data/reference/'
abadir <- '/mnt/kwiat/well/human/george/GWAS_painting/abacus/'
chroms <- 1:22

for(pop in ethnics)
{
  ## LOAD DATA
  country <- countries[ethnics == pop]
  vitlengths <- c()
  ancalleles <- c()
  for(chrom in chroms)
  {
    if(chrom < 10)
    {
      chrom <- paste0("0",chrom)
    } else
    {
      chrom <- as.character(chrom)
    }
    print(paste(pop,chrom))
    infile <- paste0(indir,pop,"GWASphasedChrom", chrom, "medianlengthsByAllele.txt.gz")
    tmp <- read.table(infile,header = T, as.is = T, stringsAsFactors = F)
    infile <- paste0(abadir,country,"/",pop,"/",chrom,"/",pop,"GWASphasedChrom", chrom, ".bayesfactor")
    abacus <- read.table(infile,header = F, as.is = T)[,3:4]
    colnames(abacus) <- c("abacus.bayesfactor","abacus.effectsize")
    vitlengths <- rbind(vitlengths,cbind(chrom,tmp,abacus))
    infile <- paste0(ancdir,"AncestralAllelesChrom", chrom, ".txt")
    tmp <- read.table(infile,header = T, as.is = T, sep = " ",stringsAsFactors = F, comment.char = "#")
    ancalleles <- rbind(ancalleles,cbind(chrom,tmp))
  }

  write.table(vitlengths,
              paste0('data/',pop,'medianviterbilengthcomps.txt'),quote = F, row.names = F, col.names = T)
}

write.table(ancalleles,paste0('data/AncestralAllelesviterbilengthcomps.txt'),quote = F, row.names = F, col.names = T)

############
# snps <- c()
# for(chrom in chroms)
# {
#   if(chrom<10) chrom <- paste0("0",chrom)
#   tmp <- read.table(paste0('/mnt/kwiat/well/human/george/GWAS_painting/abacus/Gambia/FULA/',
#                            chrom,'/FULAGWASphasedChrom',
#                            chrom,".map"), header = F, as.is = T)
#   snps <- rbind(snps,tmp)
#   
# }
############
ancalleles <- read.table(paste0('data/AncestralAllelesviterbilengthcomps.txt'),header = T)
ancalleles <- tbl_df(ancalleles)
## FIND ANC ALLELES
switches <- rep(0,nrow(ancalleles))
switches[ancalleles$ancestral_alleleA==toupper(ancalleles$alleleA)] <- 1
switches[ancalleles$ancestral_alleleA==toupper(ancalleles$alleleB)] <- 2
## PROCESS DATA

allpops <- cbind(ancalleles[,c(1,5,3,2,6,7,9)],switches)
for(pop in ethnics)
{
  print(paste('plotting',pop))
  vitlengths <- read.table(paste0('data/',pop,'medianviterbilengthcomps.txt'), header = T, as.is = T)
  colnames(vitlengths)<-gsub("median","mean",colnames(vitlengths))
  
  vitlengths <- tbl_df(vitlengths)

  ## SWITCH ALLELES TO ANC/DERIVED/UNKNOWN
  vitlengths <- mutate(vitlengths,ancswitch = switches)
  
  ## SWITCH FREQS FOR NON-DERIVED 
  newfreq <- vitlengths$a0freq
  newfreq[switches == 2] <- vitlengths$a1freq[switches == 2]
  newfreq[switches == 0] <- NA
  
  newa0 <- vitlengths$a0mean
  newa0[switches == 2] <- vitlengths$a1mean[switches == 2]
  newa0[switches == 0] <- NA
  
  newa1 <- vitlengths$a1mean
  newa1[switches == 2] <- vitlengths$a0mean[switches == 2]
  newa1[switches == 0] <- NA

  ## SWITCH FREQ AND MEANS DEPENDING ON WHICH ALLELE IS ANCESTRAL
  vitlengths <- mutate(vitlengths,
                       newa0freq = newfreq,
                       newa0mean = newa0,
                       newa1mean = newa1)
  
  ## FOR NOW LET'S DROP SNPS WHERE WE DON'T HAVE ANC/DERIVED
  ## this is about 112,000 SNPs
  vitlengths <- filter(vitlengths, !is.na(newa0mean))
  ## AND ALSO DROP MONOMORPHIC SNPS
  vitlengths <- filter(vitlengths, newa0freq != 0 & newa0freq != 1)
  
  #########
  # ## DIAGNOSTIC PLOT
  # png(paste0('figures/',pop,'lengthsAndAlleleFreqs.png'), width = 1000, height = 1000, res = 150)
  # plot((1-vitlengths$newa0freq),vitlengths$newa1mean/vitlengths$newa0mean,
  #      pch = 20, xlab = "derived allele freq",ylab = "mean derived/mean ancestral", ylim=c(0,10))
  # dev.off()
  #########
  
  ## AND SNPS WITH MAF< 0.05 | >0.95
  vitlengths <- filter(vitlengths, newa0freq >= 0.05)
  vitlengths <- filter(vitlengths, newa0freq < 0.95)

  ## ADD BINS FOR LENGTHS
  vitlengths <-  mutate(vitlengths,
                        bins1 = cut(newa0freq, seq(0,1,0.01), include.lowest = T, right = F),
                        bins5 = cut(newa0freq, seq(0.05,0.95,0.05), include.lowest = T, right = F))
  ## ADD LOG RATIOS
  vitlengths <- mutate(vitlengths,
                       logratio = log(newa0mean/newa1mean))
  
  by_bin <- group_by(vitlengths, bins5)
  bin_summary <- summarise(by_bin,
                           count = n(),
                           mean0 = mean(newa0mean, na.rm = T),
                           mean1 = mean(newa1mean, na.rm = T),
                           sd0 = sd(newa0mean, na.rm = T),
                           sd1 = sd(newa1mean, na.rm = T),
                           meanpls = mean(logratio, na.rm = T),
                           sdpls = sd(logratio, na.rm = T),
                           meanabacus = mean(abacus.effectsize, na.rm = T),
                           sdabacus = sd(abacus.effectsize, na.rm = T),
                           freqratio = mean(newa0mean/newa1mean, na.rm = T),
                           rng = IQR(newa0mean, na.rm = T))

  # ## PLOT MEAN ANCESTRAL VERSUS DERIVED LENGTHS PER BIN
  # plot(bin_summary$mean0,bin_summary$mean1,
  #      pch = 21, col = "black", bg = plot_col,
  #      xlab = "mean ancestral length per bin",
  #      ylab = "mean derived length per bin",
  #      main = "mean length comparison per bin",
  #      axes = F, xlim=c(12.8,13.4), ylim= c(12.8,13.4))
  # axis(1)
  # axis(2,las = 2)
  # abline(a=0,b=1, lwd = 2, col = "red")
  # box()

  # ## PLOT ALLELE FREQUENCY HISTOGRAMS
  # png(paste0("figures/",pop,"viterbiLengthDistibutionsByFreq100.png"), 
  #     height = 1500, width  = 1500, res = 150)
  #   p <-  ggplot(vitlengths, aes(x = newa0mean)) +
  #     geom_histogram(fill=plot_col, bins = 50 ) +
  #     facet_wrap(~bins1, scales ='free') +
  #     xlim(0,5)  +
  #     labs(title=paste0(pop," genome-wide viterbi length distibutions split by allele frequency")) 
  #   p
  # dev.off()
  # 
  # ## PLOT RAW LOG RATIOS
  # png(paste0("figures/",pop,"medianviterbiLogRatioByFreq.png"),
  #     height = 1500, width  = 1500, res = 150)
  #   p1 <-  ggplot(vitlengths, aes(x = logratio)) +
  #     geom_histogram(fill=plot_col, bins = 50 ) +
  #     facet_wrap(~bins1, scales ='free') +
  #     labs(title=paste0(pop," genome-wide log ratio of ancestral/derived lengths"))
  #   p1
  # dev.off()

  ## WE NOW HAVE A QUANTITY THAT IS SOMEWHAT EQUIVALENT TO THE
  ## raw SDS OF PRITCHARD
  ## LET'S CALL IS THE PAINTING LENGTH SCORE -- PLS
  ## WE WANT TO STANDARDISE THIS BY THE MEAN AND SD OF ALL SNPS
  ## IN THE SAME FREQUENCY BIN
  
  ## ADD BIN MEANS AND SD TO TABLE
  binsmean <- bin_summary$meanpls[match(vitlengths$bins5,bin_summary$bins5)]
  binssd <- bin_summary$sdpls[match(vitlengths$bins5,bin_summary$bins5)]
  binsabamean <- bin_summary$meanabacus[match(vitlengths$bins5,bin_summary$bins5)]
  binsabasd <- bin_summary$sdabacus[match(vitlengths$bins5,bin_summary$bins5)]
  vitlengths <- mutate(vitlengths,
                       binmean = binsmean,
                       binsd = binssd,
                       binabacusmean = binsabamean,
                       binabacussd = binsabasd)
  ## NOW COMPUTE STAN PLS AND STAN ABACUS
  vitlengths <- mutate(vitlengths,
                       stanpls = (logratio-binmean)/binsd,
                       stanabacus = (abacus.effectsize-binabacusmean/binabacussd))
                       
  
  # ## looks very normal
  # hist(vitlengths$stanpls,breaks = 100)
  # png(paste0("figures/",pop,"medianSTANabacusbayesfactorByFreq.png"),
  #     height = 1500, width  = 1500, res = 150)
  #   p1 <-  ggplot(vitlengths, aes(x = stanabacus)) +
  #     geom_histogram(fill=plot_col, bins = 50 ) +
  #     facet_wrap(~bins5, scales ='free') +
  #     labs(title=paste0(pop," genome-wide ABACUS bayes factor split by freq"))
  # 
  #   p1
  # dev.off()
   
  ## AND PVALUE?
  vitlengths <- mutate(vitlengths,
                       pvalue=2*pnorm(-abs(stanpls)))
  ## ADD TO RESULTS TABLE
  popps <- rep(NA,nrow(allpops))
  popps[match(vitlengths$snp,allpops$rsid)] <- vitlengths$pvalue
  allpops <- cbind(allpops,popps)
  colnames(allpops)[ncol(allpops)] <- paste0(pop,".p")
  popps <- rep(NA,nrow(allpops))
  popps[match(vitlengths$snp,allpops$rsid)] <- vitlengths$stanpls
  allpops <- cbind(allpops,popps)
  colnames(allpops)[ncol(allpops)] <- paste0(pop,".pls")
  
}    

### PLOT HERE

n_plots <- length(ethnics)+1
png("figures/GWASpaintingPLSmanhattenMedian.png",height = 1600, width = 2200, res = 150)
  layout(matrix(c(1:(2*(n_plots))),n_plots,2, byrow = T),
         widths = c(16,4))
  pnt_cols <- rep("darkblue",nrow(allpops))
  pnt_cols[sapply(allpops$chrom,as.numeric)%%2 == 1] <- "lightblue"
  ## PLOT MANHATTEN
  pos2 <- rep(NA,nrow(allpops))
  labs <- sapply(unique(allpops$chrom),as.numeric)
  nchr <- length(chroms)
  chromorder <- sapply(allpops$chrom, as.numeric)
  lastbase <- 0
  ticks <- NULL
  for (i in labs)
  {
    if (i == 1)
    {
      pos2[chromorder == i] <- allpops$pos[chromorder == i]
    } else
    {
      lastbase <- lastbase + tail(allpops$pos[chromorder==(i-1)],1)
      pos2[chromorder == i] <- allpops$pos[chromorder == i] + lastbase
    }
    ticks = c(ticks, pos2[chromorder==i][floor(length(pos2[chromorder==i])/2) + 1])
  }
  
  xmax <- ceiling(max(pos2) * 1.03)
  xmin <- floor(max(pos2) * -0.03)
  y_max <- 15 #max(-log10(vitlengths$pvalue))
  
  
  for(pop in ethnics)
  {
    plot_col <- pcols[ethnics == pop]
    par(mar=c(2,6,0.5,2))
    plot(pos2,rep(0,nrow(allpops)),
         ylab="",xlab="",
         xaxt="n",type="n",yaxt="n",
         ylim=c(0,y_max),xlim=c(xmin,xmax))
    mtext(2,text=expression(italic(-log[10]~P~value)), line= 2.5)
    mtext(4,las=0,text=pop,line=-2, col = pcols[ethnics == pop])
    for(i in unique(labs)) abline(v=tail(pos2[chromorder<=i],1),lty=2)
    axis(2,las=2) 
    mtext(3,text=chroms,at=ticks,  line = -2, col = "grey")
    points(pos2,-log10(allpops[,paste0(pop,".p")]),col = pnt_cols, pch = 20,cex = 0.5)
    abline(h=7,lty=2,col= "red")
    
    par(mar=c(2,4,0.5,1))
    histplot <- allpops[,paste0(pop,".pls")]
    histplot <- histplot[!is.na(histplot)]
    hist(histplot, breaks = 200,
         axes = F, xlab="", ylab= "", col = plot_col,
         main = "", xlim=c(-5,10), border = NA)
    x_at <- seq(-5,10,2.5)
    axis(1, line = -0.5, at = x_at, labels = x_at)
    axis(2,las=2, xpd = T)
    legend("topright","standardised PLS\n[painting length score]", bty = "n")
  }
  
  fishers <- -2*apply(log(allpops[,paste0(ethnics,".p")]),1,sum,na.rm = T)
  k <- length(ethnics)-apply(allpops[,paste0(ethnics,".p")],1,function(x) sum(is.na(x)))
  fishersp <- mapply(pchisq,q=fishers,df=2*k)
  fishersp <- -log10(fishersp)
  par(mar=c(2,6,0.5,2))
  plot(pos2,rep(0,nrow(allpops)),
       ylab="",xlab="",
       xaxt="n",type="n",yaxt="n",
       ylim=c(0,150),xlim=c(xmin,xmax))
  mtext(2,text=expression(italic(fishers~Chi^2)), line= 2.5)
  mtext(4,las=0,text="FISHERS META",line=-2, cex = 0.75)  
  for(i in unique(labs)) abline(v=tail(pos2[chromorder<=i],1),lty=2)
  axis(2,las=2) 
  mtext(3,text=chroms,at=ticks,  line = -2, col = "grey")
  points(pos2,fishers,col = pnt_cols, pch = 20,cex = 0.5)
  #abline(h=7,lty=2,col= "red")
  
  
dev.off()


write.table(allpops,file='/mnt/kwiat/well/human/george/GWAS_painting/length_analysis/GWASmedianPaintingLengthScores.txt',
            row.names = F, col.names = T, quote = F)

tmp <- cbind(allpops,fishers)
tmp <- tmp[tmp$chrom == 16,]
head(tmp[order(-log10(tmp$FULA.p), decreasing = T),],10)
head(tmp[order(tmp$fishers, decreasing = T),],10)



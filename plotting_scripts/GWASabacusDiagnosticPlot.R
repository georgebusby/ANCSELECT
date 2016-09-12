#######################################################
###  SCRIPT TO RUN DIAGNOSTIC PLOT WITH GWAS INDS   ###
#######################################################


## LET'S LOOK AT LCT IN FULA

hapsfile <- '/mnt/kwiat/well/human/george/GWAS_painting/abacus/Gambia/FULA/02/FULAGWASphasedChrom02.hap'
mapfile <- '/mnt/kwiat/well/human/george/GWAS_painting/abacus/Gambia/FULA/02/FULAGWASphasedChrom02.map'
vitfile <- '/mnt/kwiat/well/human/george/GWAS_painting/abacus/Gambia/FULA/02/FULAGWASphasedChrom02'

## LOAD FILES
haps <- read.table(hapsfile,header = F)
map <- scan(mapfile)
snps <- cbind(haps[,1:5],map)
colnames(snps) <- c("SNPID","rsid","pos","a0","a1","map")
haps <- haps[,6:ncol(haps)]

## SELECT SNPS
snpname <- "rs13006092"
snprange <- c(132e6,142e6)
snpsindex <- which(snps$pos>=snprange[1]&snps$pos<=snprange[2])

## SELECT HAPS TO PLOT
hap_plot <- haps[snpsindex,]
focalsnp <- which(snpsindex==which(snps$rsid==snpname))
anchaps <- hap_plot[focalsnp,] == 0

## SEPARATE HAPS AND KEEP ORDER
## ANCESTRAL HAPS
hap_plot0 <- hap_plot[,anchaps]
## ORDER
neword <- hclust(dist(t(hap_plot0)))$order #[(focalsnp-50):(focalsnp+50),]
hap_plot0 <- hap_plot0[,neword]
ancord <- which(anchaps)[neword]

## ADD A GAP
gap <- 50

## DERIVED HAPS
hap_plot1 <- hap_plot[,!anchaps]
neword <- hclust(dist(t(hap_plot1)))$order #[(focalsnp-50):(focalsnp+50),]
hap_plot1 <- hap_plot1[,neword]
derord <- which(!anchaps)[neword]

## COMBINE TOGETHER
hap_plot01 <- cbind(as.matrix(hap_plot0),
                    matrix(0.5,nr=nrow(hap_plot0),nc = gap),
                    as.matrix(hap_plot1))

# ## LOAD LENGTHS AND COPIER FOR HAPS
library(data.table)
firstrow <- snpsindex[1]
numrows <- 1+(snpsindex[length(snpsindex)]-snpsindex[1])

for(i in 0:(ncol(haps)-1))
{
  print(i)
  vit <- fread(paste0(vitfile,"_",i,".viterbi"), skip = firstrow - 1, nrows = numrows)
  if(i == 0)
  {
    vitlengths <- vit$V4
    vitcopy <- vit$V3
  } else
  {
    vitlengths <- cbind(vitlengths,vit$V4)
    vitcopy <- cbind(vitcopy,vit$V3)
  }
}


## PLOT
png()

## PLOT UNDERLYING HAPS
pos_plot <- snps$pos[snpsindex]
image(pos_plot,
      1:ncol(hap_plot01),
      hap_plot01,
      col = c("white","white","grey"),
      axes = F, xlab = "", ylab = "")

## REORDER VITCOPY TO BE SAME AS PLOT_HAP
realord <- c(ancord,derord)
vitcopy2 <- vitcopy3 <- vitcopy[,realord]

## NOW WE WANT TO ADD THE VITERBI PATHS ON TOP
nocolour <- rgb(0,0,0,0)
vithaps <- matrix(0,nr=dim(vitcopy2)[1],nc=dim(vitcopy2)[2])

## NEED TO COLOUR LENGTHS AROUND FOCAL SNP
## CONVERT LENGTH AT SNP TO SNP DISTANCE
tmp <- apply(vitcopy2,2,function(x){x==x[focalsnp]})
vithaps[,1:sum(anchaps)][tmp[,1:sum(anchaps)]] <- 1
vithaps[,sum(anchaps):ncol(vithaps)][tmp[,sum(anchaps):ncol(vithaps)]] <- 2

vithaps <- cbind(vithaps[,1:sum(anchaps)],
                 matrix(0,nr=nrow(hap_plot0),nc = gap),
                 vithaps[,(sum(anchaps)+1):ncol(vitcopy2)])

image(pos_plot,
      1:ncol(vithaps),
      vithaps,
      col = c(nocolour, "blue", "red"),
      axes = F, xlab = "", ylab = "",
      add = T)
xat <- pretty(pos_plot)
axis(1,at = xat, labels = xat/1e6, xpd = T)

## NOW PLOT BY ALLELE THAT COPIER IS CARRYING

hap_plot2 <- hap_plot[,realord]
sameallele <- hap_plot2[focalsnp,vitcopy3[focalsnp,]] == hap_plot2[focalsnp,]
tmp <- apply(vitcopy3,2,function(x){x==x[focalsnp]})

vithaps <- matrix(0,nr=dim(vitcopy2)[1],nc=dim(vitcopy2)[2])
wrongallele <- sameallele[1:sum(anchaps)]
vithaps[,!wrongallele][tmp[,!wrongallele]] <- 1
vithaps[,wrongallele][tmp[,wrongallele]] <-2
newvithaps <- cbind(vithaps[,!wrongallele[order(vitlengths[focalsnp,!wrongallele])]],
                    vithaps[,wrongallele[order(vitlengths[focalsnp,wrongallele])]])

wrongallele <- sameallele[sum(anchaps):ncol(vithaps)]
vithaps[,!wrongallele][tmp[,!wrongallele]] <- 1
vithaps[,wrongallele][tmp[,wrongallele]] <- 2
newvithaps <- cbind(newvithaps,
                    vithaps[,!wrongallele],
                    vithaps[,wrongallele])
vithaps <- cbind(newvithaps[,1:sum(anchaps)],
                 matrix(0,nr=nrow(hap_plot0),nc = gap),
                 newvithaps[,(sum(anchaps)+1):ncol(vitcopy2)])

## PLOT UNDERLYTING HAPS
image(pos_plot,
      1:ncol(hap_plot01),
      hap_plot01,
      col = c("white","white","grey"),
      axes = F, xlab = "", ylab = "",
      main = paste("red = same SNP; blue = different SNP"))

image(pos_plot,
      1:ncol(vithaps),
      vithaps,
      col = c(nocolour, "blue", "red"),
      axes = F, xlab = "", ylab = "",
      add = T)
xat <- pretty(pos_plot)
axis(1,at = xat, labels = xat/1e6, xpd = T)









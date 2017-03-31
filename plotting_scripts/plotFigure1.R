### SCRIPT TO CALL THE VARIOUS PLOTTING SCRIPTS TO PRODUCE A SINGLE PLOT ###
setwd("~/repos/ANCSELECT/")

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
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
#popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
#pca_file <- "data/Africa300KPCS.txt"
#tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
#mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)
TAB <- paste(rep(" ", 4), collapse = "")

pdf("figures/LocalAncestryOverviewFig1.pdf",width=10,height=5)
    par(mar=c(0.5,0.5,0.5,0.5))
    layout(matrix(c(1,3,
                    1,4,
                    2,5),3,2,byrow=T),
           heights=c(1,2,2),widths=c(3,7))

    panel_desc <- 3.45 ## position of the panel descriptions in the right margin
    ###########################################################
    ## 00 MAP AND LEGENDS
    transparent_cols <- FALSE ## use transparent colours for malariagen countries
    plot_equator <- FALSE ## should equator and tropics be plotted
    plot_ancestry_regions <- TRUE ## should countries be coloured by ancestry region or individually?
    map_xlim<-c(-20,50) ## X-AXIS LIMITS OF MAP
    map_ylim<-c(-32,28) ## Y-AXIS LIMITS OF MAP
    #regions <- unique(popkey$RegionM) ## regions
    pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
    ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                                "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                                "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
    map_ocean_col <- "white" ## colour of ocean in map
    map_border <- "grey10"
    map_miss_col <- "grey" ## colour of non-focal countries in map
    pt_cex <- 2
    pt_lwd <- 0.75
    plot_legends <- FALSE
    plot_letter <- FALSE
    plot_points <- FALSE
    plot_box <- F
    source("~/repos/admixture_in_africa/plotting_scripts/africamap.R")
    ## ADD POINTS FOR COUNTRY SAMPLING POSITIONS
    latlongs <- read.csv(latlong_file,header=T)
    points(latlongs$Long,latlongs$Lat,pch=21,cex=0.5,
           col = "grey50", bg = "slateblue3")    
    mtext(3,line=-2,text=" a",adj=0, cex = 2)
    ## ADD A LEGEND
    plot(0,0,type = "n", axes = F, xlab = "", ylab = "")
    legend("topleft",legend=c(gsub("Africa ","",gsub("\\_"," ",regions)),"populations from \nBusby et al 2016"),
           pch=c(rep(22,length(regions)),21),
           pt.bg=c(pcolshex,"slateblue3"),
           col=c(pcolshex,"grey50"),
           bty="n",title="Ancestry Regions", 
           pt.cex = c(rep(2,length(regions)),1),
           pt.lwd = 1,
           ncol = 2)

    #########################################################
    ## 01 PLOT A PAINTED CHROMOSOME
    library("h5")
    ## FUNCTIONS
    hap2sampleindex <- function(hap,nsamps=10){
        ## finds the first sample index for a haplotype
        sample <- (hap*nsamps)-(nsamps-1)
        return(sample)
    }

    datafile <- h5file('/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5')
    ## CHROMOSOME 2
    chrom <- "22"
    ## GET MAP AND POSITION INFO
    map <- data.frame(readDataSet(datafile[paste0("/paintings/chrom",chrom,"/map")]))
    colnames(map) <- c("position","recrate")
    ## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
    psamples <- readDataSet(datafile[paste0("/paintings/samples/individuals")])
    colnames(psamples) <- c("ind","region","X")
    
    ## PLOT TEN DIFFERENT INDS
    inds2plot <- c(paste0("FULA_S", 1:10))
    psamplesind <- (1:nrow(psamples))[psamples[,"ind"] %in% inds2plot ]
    
    ## 2 haps per sample!!
    psampleshap <- hap2sampleindex(psamplesind,2)
    psamplesindsamp <- hap2sampleindex(psampleshap)
    
    analysis <- "nonlocal"
    #tmp <- psamplesindsamp:(psamplesindsamp+9) ## FOR ONE INDIVIDUAL
    tmp <- psamplesindsamp
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
        
    ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
    paintedchromregprop <- matrix(0,nc=length(ancreg_list),nr=nrow(paintedchromreg))
    for(i in 1:nrow(paintedchromreg))
    {
      #print(i)
      tmp <- table(paintedchromreg[i,])
      paintedchromregprop[i,as.numeric(names(tmp))] <- tmp
    }
    
    ###########################################################
    ## PLOT ACTUAL HAPLOTYPES
    chromlength <- as.numeric(as.character(map$recrate))
    chrompos <- as.numeric(as.character(map$position))
    chromposI <- c(diff(chrompos),0)
        
    chromplot <- paintedchromreg
    chromcolsbreaks <- sort(unique(unlist(apply(paintedchromreg,2,unique))))
    chromcols <- pcolshex[chromcolsbreaks]
    par(mar=c(0,1,2,7))
    image(chrompos,
          1:ncol(chromplot),
          (chromplot),
          col=chromcols,breaks = c(0,chromcolsbreaks),
          axes=F,xlab="",ylab="")
    mtext(3,line=0,text=" b",adj=0, cex = 2)
    ## PLOT AXIS
    # xatlab <- pretty(chrompos)
    # #xat <- seq(0,sum(chrompos), length.out = length(xatlab))
    # #xat[length(xat)] <- sum(chromposI)
    # axis(1,at=xatlab,labels=xatlab/1e6)
    # mtext(1,text=paste("position on chromosome",as.numeric(chrom), "(Mb)"),line=2,cex=1)
    mtext(4, text = "10 randomly\nselected painted\nFula haplotypes",
          las = 2, adj = 0.5, line = panel_desc, cex = 0.75)

    ####################################################################################
    ## NOW SHOW PAINTINGS ACROSS ALL FULAI
    psamplesind <- (1:nrow(psamples))[psamples[,"region"] == "FULAI"]
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
    ## SUBSAMPLE SNPS
    snpsample <- round(seq(1,nrow(chromplot), length = 1000))
    ####################################################################################
    ## ALL PAINTED CHROMSOMES    
    par(mar=c(0,1,2,7))
    image(chrompos[snpsample],
          1:ncol(chromplot),
          (chromplot[snpsample,]),
          col=chromcols,breaks = c(0,chromcolsbreaks),
          axes=F,xlab="",ylab="")
    # mtext(3,
    #       text=paste("144 Fulani chromosomes painted by different donor region"), #,"[",nrow(chromplot), "SNPs on chromosome", chrom,"]"
    #       line=0.5,cex=1.5, adj=0.5)
    xatlab <- pretty(chrompos)
    xat <- seq(min(chrompos),max(chrompos), length.out = length(xatlab))
    # #xat[length(xat)] <- sum(chromposI)
    # axis(1,at=xat,labels=xatlab/1e6)
    
    mtext(3,line=0,text=" c",adj=0, cex = 2)
    mtext(4, text = "All painted\nFula haplotypes",
          las = 2, adj = 0.5,line = panel_desc, cex = 0.75)
    
    
    ####################################################################################
    ## PROPORTIONS    
    chromplot <- paintedchromregprop[,8:1]
    par(mar=c(3,1,2,7))
    barplot(t(chromplot),
            width=chromposI,
            col=rev(pcolshex),xaxs="i",yaxs="i",
            space=0,axes=F,xaxt="n",border=NA,
            xlim=c(0,sum(chromposI)),horiz=F,
            xlab="")
    xatlab <- pretty(chrompos)
    xat <- seq(0,sum(chromposI), length.out = length(xatlab))
    #xat[length(xat)] <- sum(chromposI)
    axis(1,at=xat,labels=xatlab/1e6)
    mtext(1,text=paste("position on chromosome",as.numeric(chrom), "(Mb)"),line=2,cex=1)
    mtext(3,line=0,text=" d",adj=0, cex = 2)
    mtext(4, text = "Ancestry\nproportions\nacross all\npainted \nFula haplotypes",
          las = 2, adj = 0.5, line = panel_desc, cex = 0.75)
 
    ## ADD RECT HIGHLIGHTING A DEVIATION
    xleft <- (21.75/40)*sum(chromposI)
    xright <- (23.25/40)*sum(chromposI)
    rect(xleft,0,xright,1250, border = "black", lwd = 1)    
    midpoint  <- mean(c(xleft,xright))
    arrows(midpoint,1450,midpoint,1250, length = 0.05, lwd = 1, angle = 20, xpd = T)    
    mtext(3,text = "Is Eurasian ancestry (dark yellow) significantly deviated from expectations?",
          cex = 0.5)
         
dev.off()

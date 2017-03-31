### SCRIPT TO CALL THE VARIOUS PLOTTING SCRIPTS TO PRODUCE A SINGLE PLOT ###
setwd("~/repos/ANCSELECT/")
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
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
#popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
#pca_file <- "data/Africa300KPCS.txt"
#tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
#mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")
datafile <- h5file('/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5')


## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$RegionM <- gsub("Afro-Asiatic","Afroasiatic",popkey$RegionM)
TAB <- paste(rep(" ", 4), collapse = "")


#####
## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
  ## finds the first sample index for a haplotype
  sample <- (hap*nsamps)-(nsamps-1)
  return(sample)
}



#####


pdf("figures/FulaLocalAncestryDeviationsFig2malaria.pdf",width=10,height=5)
    layout(matrix(c(4,9,
                    3,8,
                    1,6,
                    2,7,
                    5,10),5,2,byrow=T),
           heights=rep(c(1,0.5,1,1,1.5),1),
           widths=c(1,1))

    ##############################################
    ## SET UP
    pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
    ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                                "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                                "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

    analysis <- "nonlocal"
    
    for(run in 1:2)
    {
    
      if(run == 1)
      {
        chrom <- "02"
        plot_range <- c(120e6,130e6)
        pop <- "FULAI"
        #pop <- "MKK"
        panel_letter <- "a"
        gene_highlight <- c("CNTNAP5")
      }

      if(run == 2)
      {
        chrom <- "04"
        plot_range <- c(142e6,148e6)
        pop <- "FULAI"
        #pop <- "AMHARA"
        panel_letter <- "b"
        gene_highlight <- c("GYPE","GYPA","GYPB")
      
      }

      #######################################################
      ### PLOT LOGISTIC AND MVN
      library(dplyr) ## for binding dataframes with difference numbers of columns; bind_rows
      #in_dir <- '/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/'
      in_dir <- '/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/'
      pops <- read.table('/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/analysislists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP.idfile.txt')
      pops <- as.character(levels(pops[,2]))
     
      chrom2 <- as.numeric(chrom)
      in_file <- paste0(in_dir,pop,'nolocalChrom',chrom,'.ancestryselectionIII.gz')
      in_file2 <- paste0(in_dir,pop,'nolocalChrom',chrom,'.ancestryselectionIIIa.gz')
      liks <- read.csv(in_file,header=T,as.is=T)
      mvn <- read.csv(in_file2,header=T,as.is=T)
      pcols <- grep(".RC.P",colnames(liks))
      plot_region <- liks$pos>=plot_range[1]&liks$pos<=plot_range[2]
      liks <- liks[plot_region,]
      ymax <- max(liks[,pcols])    
#      ylims <- c(0,ymax)
 
      ylims <- c(0, 15)    
      ## FIND THE RANGE OVER WHICH THE PVAL IS HIGH
      maxp <- which.max(apply(liks[,pcols],2,max))
      maxp <- c(which.max(liks[,pcols[maxp]]), maxp)
      pvalue <- liks[maxp[1],pcols[maxp[2]]]
      minp <- 1
      
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
      
      gene_region <- c(liks$pos[pos1],liks$pos[pos2])
      
      par(mar=c(3,6,1,1))
      plot(0,0,
           xlim=range(plot_range), type = "n", ylim = ylims,
           axes = F, ylab = expression(-log[10]~italic(P)), xlab = "",
           xaxs = "i", yaxs = "i")
      abline(v=gene_region, col = "grey10", lty = 1, xpd = T)
      
      ## ADD MVN LINE
      points(mvn$pos,mvn$MVNp, col = "black", type = "S", lwd = 2)
      
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
  
      #######################################################
      ### PLOT BETA COEFFS
      pcols <- grep(".RC.beta",colnames(liks))
      #ylims <- range(pretty(range(liks[,pcols])))
      ylims <- c(-4.5,2.5)
      
      par(mar=c(0,6,0,1))
      plot(0,0,
           xlim=range(plot_range), type = "n", ylim = ylims,
           axes = F, ylab = expression(beta~Anc.~deviation), xlab = "",
           xaxs = "i", yaxs = "i")
      abline(h=0, lty = 2)
      abline(v=gene_region, col = "grey10", lty = 1, xpd = T)
      
      
      for(i in pcols)
      {
        reg <- gsub("\\.","-",gsub("\\.RC.beta","",colnames(liks)[i]))
        linecol <- pcolshex[regions == reg]
        points(liks$pos,liks[,i], col = linecol, type = "S", lwd = 2)
      }
      yat <- pretty(range(ylims), n= 4)
      axis(2, at = yat, labels = yat, las = 2, xpd = T)
      
      #########################################################
      ## 01 PLOT A PAINTED CHROMOSOME
      
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
      # test_rev_region <- chrompos>=gene_region[1]&chrompos<=gene_region[2]
      # tmp <- paintedchrom[test_rev_region,]
      # 
      ####################################################################################
      ## SUBSAMPLE SNPS
      #snpsample <- round(seq(1,nrow(chromplot), length = 1000))
      ####################################################################################
      ## PROPORTIONS    
      chromplot <- paintedchromregprop[,8:1]/rowSums(paintedchromregprop)
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
      par(new = T)
      plot(0,0,axes = F, ylim=c(0,1), type = "n",
           xlab = "", ylab = "", xlim=plot_range, xaxs = "i", yaxs = "i")
      abline(v=gene_region, col = "grey10", lty = 1, xpd = T)
      
           
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
      
      par(mar=c(0,6,1,1))
      plot(gm$position[gm$chromosome==chrom],gm$COMBINED_LD_rate[gm$chromosome==chrom]*1e6,
           xlim=range(plot_range), type = "S", ylim = c(0,100),
           axes = F, ylab = "Rec Rate\n(cM/Mb)", xlab = "", xaxs = "i")
      axis(2, las = 2)
      points(gm$position[gm$chromosome==chrom],gm$YRI_LD_rate[gm$chromosome==chrom]*1e6,
             type = "S", col = "blue", lty = 2)
      points(gm$position[gm$chromosome==chrom],gm$African_Enriched_rate[gm$chromosome==chrom]*1e6,
             type = "S", col = "red", lty = 3)
      legend("topleft",legend=c("HAPMAP","YRI", "Africa enriched"),
             col = c("black","blue","red"), lty = c(1,2,3), bty = "n", horiz = T)
  
      mtext(2,at=100,text=panel_letter,adj=0, cex = 1.5, las = 2, line = 5)          
      #######################################################
      ### PLOT GENES
      source("~/repos/glycophorins/external_software/plot_genes.R")
      genes = load.genes( "/mnt/kwiat/data/1/galton/malariagen/human/reference/genome-mysql.cse.ucsc.edu/2014-03-18/UCSC_hg19_2014-03-18_refGene.tsv" )
      genes <- genes[genes$cdsStartStat!="unk",]
      par(mar=c(4,6,1,1))
      
      if(run == 1)
      {
        plot.genes( chromosome = chrom, region = range(pretty(gene_region)), genes, xaxt = "n" , plot.ylab = "Genes", highlight.gene = gene_highlight, switch.gene.side = TRUE)
      }
      
      if(run == 2)
      {
        plot.genes( chromosome = chrom, region = range(pretty(gene_region)), genes, xaxt = "n" , plot.ylab = "Genes", highlight.gene = gene_highlight)
      }
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

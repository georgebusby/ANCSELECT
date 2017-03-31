#######################################################
### SCRIPT TO PLOT LRT ACROSS GWAS ETHNIC GROUPS    ###
#######################################################
library("dplyr")
source("plotting_scripts/gavinManhattenFunctions.R")

#############################################
## GET DATA
snps <- read.table("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330Kphased.legend", header=F,as.is=T)
colnames(snps) <- c("chrom","rsid","pos","a0","a1")
popkey_file <- "data/PopulationKey.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$AncestryRegion <- gsub("Afro-Asiatic","Afroasiatic",popkey$AncestryRegion)
colour_table <- read.table("data/RegionColours.txt", header = T, stringsAsFactors = F, comment.char = "")
ancreg_list <- colour_table$AncestryRegion
pcolshex <- colour_table$Colour
#############################################
## GET LRT RESULTS
pop <- "FULAI"
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                            "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afroasiatic",
                            "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

##poplist <- popkey$Ethnic_Group
poplist <- c("JOLA","FULAI","JUHOAN")
## start at OROMO
for(pop in poplist)
{
   
    y_max <- 20 #round(max(fsts),2)
    maxx <- 18
    main.signifp <- 8
  ## MKK?? AN ALL AFOR ASIATIC/NILO SAHARANS
    print(paste("plotting",pop))
    if(pop == "SEMIBANTU") pop = "SEMI.BANTU"
    if(pop == "HERERO") pop = "SWBANTU"
    if(pop == "JUHOANSI") pop = "JUHOAN"
    if(pop == "GUIGANA") pop = "GUIGHANAKGAL"
    if(pop == "TIGRAY") pop = "TYGRAY"
    chroms <- 1:22
    # in_root <- paste0("/mnt/kwiat/data/2/bayes/users/george/popgen/analysis3/ancestry_selection/",pop,"nolocalChrom")
    in_root <- paste0("/mnt/kwiat/well/human/george/copy_selection2/copy_selection/output/",pop,"nolocalChrom")
    in_suff <- ".ancestryselectionIII.gz"
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
    if(pop == "SEMI.BANTU") pop = "SEMIBANTU"
    if(pop == "SWBANTU") pop = "HERERO"
    if(pop == "JUHOAN") pop = "JUHOANSI"
    if(pop == "GUIGHANAKGAL") pop = "GUIGANA"
    if(pop == "TYGRAY") pop = "TIGRAY"
    
    ## PRELIMINARY PLOTS
    signifp <- main.signifp
    
    self_reg <- gsub("\\-","\\.",popkey$AncestryRegion[popkey$Ethnic_Group==pop])
    pos2 <- rep(NA,nrow(snps))
    labs <- unique(snps$chrom)
    nchr <- length(chroms)
    lastbase <- 0
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
    }
    
    xmax <- ceiling(max(pos2) * 1.03)
    xmin <- floor(max(pos2) * -0.03)
    pos2 <- rev(pos2)
    ticks <- NULL
    for(i in labs)
    {
      ticks <- c(ticks,median(pos2[pl$chrom==i]))
    }
    ################################
    ## POINT COLOURS
    cols <- rep("grey10",nrow(pl));
    cols[pl$chrom %in% seq(2,22,by=2) ] <- "grey30";
    ################################
    ## LET'S PLOT ALL ANCESTRIES ON THE SAME PLOT
    ## PLOT EVERTHING FIRST AND THEN ADD HITS
    #png(paste0("figures/",pop,"LRTmanhattanIII.png"),height = 2000, width = 2200,res = 300,pointsize = 10)
    pdf(paste0("figures/",pop,"LRTmanhattanIII.pdf"), height = 10, width = 7)
    
    layout(matrix(c(1,3,2),1,3),widths =c(3,1.5,7.5))
    
    topmar <- 3
    botmar <- 4
    par(mar=c(botmar,4,topmar,0))
    plot(rep(0,nrow(snps)),pos2,
         ylab="",xlab="",yaxs = "i",xaxs = "i",
         xaxt="n",type="n",yaxt="n",
         xlim=c(0,y_max+2),ylim=c(xmin,xmax), axes = F)
    mtext(1,text=expression(-log10(P)), las = 1, line= 2.5, cex = 0.75)
    #mtext(2,las=0,text=pop,line=-2, col = "black")
    for(i in unique(labs)) abline(h=tail(pos2[snps$chrom<=i],1),lty=2, xpd = F, col = "grey")
    axis(1,las=1)
    mtext(2,text=chroms,at=ticks,  line = 1, col = "black", las = 2, cex = 0.75)
    axis(2,at = ticks, labels = NA)
    abline(v= 0, xpd = F, lwd = 2)
    
    ## ADD GENOME-WIDE SIGNIFICANCE
    signifp <- main.signifp
    abline(v= signifp, col = "red", lty = 3, xpd = F)
    hits <- c()
    plot_allpnts <- T
    while(is.null(hits))
    {
      self_reg2 <- self_reg
      if(self_reg %in% c("East_Africa_Nilo.Saharan","East_Africa_Afroasiatic"))
      {
        self_reg2 <- c("East_Africa_Nilo.Saharan","East_Africa_Afroasiatic") 
      }
      for(reg in ancreg_list[!ancreg_list%in%gsub("\\.","-",self_reg2)])
      {
        regc <- gsub("\\-","\\.",paste0(reg,".RC.P"))
        regb <- gsub("\\-","\\.",paste0(reg,".RC.beta"))
        numdon <- gsub("\\-","\\.",paste0(reg,".num.dons"))
        prop <- gsub("\\-","\\.",paste0(reg,".prop"))
        ps <- pl[,regc]
        tmp <- pl[which(ps>signifp),c(colnames(pl)[1:6],regb,regc,numdon,prop)]
        pntcex <- rep(0.5,length(ps))
        pntcex[which(ps>=signifp)] <- 0
        if(plot_allpnts) points(ps,pos2, cex = pntcex, col = cols, pch = 19)  
        if(nrow(tmp)>0)
        {
          ## BIN HITS INTO 2.5Mb REGIONS
          tmp$pos2 <- cut(tmp$pos, seq(0,max(pl$pos),by=2.5e6))
          hitpostab <- unique(tmp[,c("chrom","pos2")])
          for(j in 1:nrow(hitpostab))
          {
            tmp2 <- tmp[tmp$chrom==hitpostab$chrom[j]&tmp$pos2==hitpostab$pos2[j],]
            tophit <- tmp2[which.max(tmp2[,regc]),]
            tmpout <- unlist(c(tophit[,c(1:3,7,8,9,10)],reg))
            hits <- rbind(hits,tmpout)
            colnames(hits)[4:8] <- c("beta","P","num.dons","prop","region")
          }
        }
      }
      signifp <- signifp - 1
    }
    
    ## unique hits
    hits <- data.frame(hits,stringsAsFactors = F)
    whichna <- which(!is.na(as.numeric(hits$chrom)))
    hits <- hits[whichna,]
    hits$chrom <- as.numeric(hits$chrom)
    hits$pos <- as.numeric(hits$pos)
    hits$P <- as.numeric(hits$P)
    hits$beta <- as.numeric(hits$beta)
    hits <- hits[order(hits$chrom,hits$pos),]
    
    ## FIND GENOMIC REGION AROUND HIT
    hits$lowerpos <- hits$upperpos <- 0
    for(j in 1:nrow(hits))
    {
     ## FIND SNPS AROUND HIT WHERE PVALUE IS >=minp
      topp <- which(pl$rsid==hits$rsid[j])
      minp <- 4
      regc <- gsub("\\-","\\.",paste0(hits$region[j],".RC.P"))
      i <- pl[topp,regc]
      posl <- topp
      tchrom <- pl$chrom[posl]
      ichrom <- pl$chrom[posu-1]
      if(length(ichrom) == 0 || is.na(ichrom)) ichrom <- 0
      while(i > minp & tchrom == ichrom)
      {
        posl <- posl - 1
        i <- pl[posl,regc]
        ichrom <- pl$chrom[posl-2]
        if(length(ichrom) == 0 || is.na(ichrom)) ichrom <- 0
        if(is.na(i)) i <- 0
      }
      i <- pl[topp,regc]
      posu <- topp
      tchrom <- pl$chrom[posu]
      ichrom <- pl$chrom[posu+1]
      if(length(ichrom) == 0 || is.na(ichrom)) ichrom <- 23
      while(i > minp & tchrom == ichrom)
      {
        posu <- posu + 1
        i <- pl[posu,regc]
        ichrom <- pl$chrom[posu+1]
        if(length(ichrom) == 0 || is.na(ichrom)) ichrom <- 23
        if(is.na(i)) i <- 0
      }
      hits$lowerpos[j] <- pl$pos[posl]
      hits$upperpos[j] <- pl$pos[posu]
    }
    ## DO ANOTHER ROUND OF PRUNING OF SIGNALS
    anydups <- which(duplicated(hits[,c("upperpos","lowerpos")]))
    if(length(anydups)>0)
    {
      dups2go <- c()
      for(i in anydups)
      {
        indices <- (i-1):i
        tmp1 <-hits[indices,]
        togo <- indices[which.min(tmp1$P)]
        dups2go <- c(dups2go,togo)
      }
      hits <- hits[-dups2go,]
    }
    
    ###########################################################################
    ## TO GUARD AGAINST SITUATIONS WHERE COPYING PROPORTIONS MIGHT BE INFLUENCED
    ## BY AN ADMIXED DONOR INDIVIDUAL, WE REMOVE HITS WHERE THE NUMBER OF DONORS
    ## CONTRIBUTING TO THE SIGNAL IS SMALL (<=3) OR WHERE A DISPROPORTIONATELY
    ## SMALL NUMBER OF INDIVIDUALS CONTRIBUTE TO THE REGION'S ANCESTRY PROPORTION
    ## WE COMPUTE THIS BY:
    ## number of unique donor individuals/total number of donors
    
    minimum.num.dons <- 3
    minimum.don.prop <- 0.1
    
    test <- which((as.numeric(hits$num.dons)>minimum.num.dons) & (as.numeric(hits$num.dons)/as.numeric(hits$prop) > minimum.don.prop))
    ## NOW WE WANT TO FIND THE CLOSEST GENES
    source("~/R/Copy/Rprojects/AfricaPOPGEN/functions/hitplots.R")
    genes <- load.genes()
    genes <- genes[genes$cdsStartStat!="unk",]
    ginfo <- matrix(0,ncol=22,nrow=0)
    colnames(ginfo) <- c(colnames(genes),"rsid","reg","P","upper","lower")
    if(nrow(hits) > 0)
    {
      for(i in 1:nrow(hits))
      {
        genes2 <- genes[ which( genes$chrom == paste0("chr",hits$chrom[i]) & genes$txEnd >= (hits$lowerpos[i]) & genes$txStart <= (hits$upperpos[i])),, drop = FALSE ]
        if(nrow(genes2)>0)
        {
          genes2 <- cbind(genes2,hits$rsid[i],hits$region[i],hits$P[i],hits$lowerpos[i],hits$upperpos[i])
          colnames(genes2) <- colnames(ginfo)
        } else
        {
          genes2 <- data.frame(t(c(rep(NA,17),hits$rsid[i],hits$region[i],hits$P[i],hits$lowerpos[i],hits$upperpos[i])),stringsAsFactors = F)
          if(!is.null(colnames(ginfo))) colnames(genes2) <- colnames(ginfo)
        }
        ginfo <- rbind(ginfo,genes2)
      }
    
      ginfo <- data.frame(ginfo, stringsAsFactors = F)
      
      ############################################################
      ## STORE AND TRIM HITS
      hit_tab <- paste0("data/",pop,"LRTmanhattanIII_hits.txt")
      ginf_tab <- paste0("data/",pop,"LRTmanhattanIII_ginfo.txt")
      
      keep <- rep(0, nrow(hits))
      keep[test] <- 1
      
      write.table(cbind(hits,keep), file=hit_tab, sep = "\t", col.names = T, row.names = F, quote = F)
      write.table(ginfo, file=ginf_tab, sep = "\t", col.names = T, row.names = F, quote = F)
      
      ###########################################################################
      ## TRIM
      hits <- hits[test,]
      ###########################################################################
      ## LIMIT TO TOP 30 HITS MAX
      n_tophits <- 25
      if(nrow(hits) > n_tophits )
      {
        hits <- hits[order(hits$P, decreasing = T)[1:n_tophits],]
        hits <- hits[order(hits$chrom,hits$pos),]
      }

      hit_rsid <- unique(as.character(hits$rsid))
      ginfo <- ginfo[which(ginfo$rsid%in%hit_rsid),]
      ############################################################
      ## FIND LENGTH OF GENE NAMES IN CHARACTGERS
      num_gene_rows <- c()
      for(i in 1:nrow(hits))
      {
        #if(!is.na(ginfo$rsid))
        gene_info <- paste(ginfo[ginfo$rsid==hits$rsid[i],c("name2")],collapse = " ")
        if(nchar(gene_info)<= 35) num_gene_rows <- c(num_gene_rows,2)
        if(nchar(gene_info)> 35 & nchar(gene_info) <=70) num_gene_rows <- c(num_gene_rows,2)
        if(nchar(gene_info)> 70 & nchar(gene_info) <=105) num_gene_rows <- c(num_gene_rows,3)
        if(nchar(gene_info)> 105& nchar(gene_info) <=140) num_gene_rows <- c(num_gene_rows,4)
        if(nchar(gene_info)> 140) num_gene_rows <- c(num_gene_rows,5)
      }
      num_gene_rows <- c(2,num_gene_rows)
      
      ## FIND SPLIT POSITIONS
      numhits <- nrow(hits)
      ## split the y axis into numhits bins
      splits <- seq(xmax,xmin,length = sum(num_gene_rows))
      splits2 <- c()
      for(i in 1:length(num_gene_rows))
      {
        index <- sum(num_gene_rows[1:i])
        splits2 <- c(splits2, splits[index])
      }
      splits <- c(splits2)
      
      midsplits <- c()
      for(i in 2:length(splits)) midsplits <- c(midsplits,mean(c(splits[(i-1)],splits[i])))
      #for(i in midsplits) abline(h=i, xpd = F, col = "red", lty = 2)
      
      ## NOW DRAW LINES FROM THE HIT TO THE EDGE AND PLOT POINTS
      if(nrow(hits) > 0)
      {
        for(j in 1:nrow(hits))
        {
          ## DRAWLINES
          topp <- which(pl$rsid==hits$rsid[j])
          x0 <- hits$P[j]
          y0 <- pos2[topp]
          x1 <- maxx
          y1 <- y0
          x2 <- y_max+2
          y2 <- midsplits[j]
          segments(x0,y0,x1,y1)
          segments(x1,y1,x2,y2)  
          posu <- hits$upperpos[j]
          posl <- hits$lowerpos[j]
          hitchrom <- hits$chrom[j]
          hits2plot <- which(pl$pos==posl&pl$chrom==hitchrom):which(pl$pos==posu&pl$chrom==hitchrom)
          hitscol <- pcolshex[ancreg_list==gsub("\\.","\\-",hits$region[j])]
          regc <- gsub("\\-","\\.",paste0(hits$region[j],".RC.P"))
          points(pl[hits2plot,regc],pos2[hits2plot],col = hitscol,pch = 18, cex = 1)
        }
      }
    }
    ###########################################################################
    ## NOW ADD A SECOND PLOT WITH ALL THE INFO ABOUT THE TOPHITS
    par(mar=c(botmar,0.5,topmar,2))
    xlims <- c(0,5)
    textpos <- c(0,0,0.65,2,1)
    texttitles <- c("Position", "Variant","P", "Gene(s) of Interest", "Genomic Region\nsize(Mb)")
    plot(rep(0,nrow(snps)),pos2,
         ylab="",xlab="",yaxs = "i",xaxs = "i",
         xaxt="n",type="n",yaxt="n",axes = F,
         xlim=xlims,ylim=c(xmin,xmax))
    
    ## ADD SPLIT LINES
    for(i in splits) abline(h=i, xpd = F)
    abline(h=splits[1], xpd = F, lwd = 2)
    
    annotext <- 0.8
    
    if(nrow(hits) > 0)
    {
      for(i in 1:nrow(hits))
      {
        text1 <- paste0(hits$chrom[i],":",round(hits$pos[i]/1e6,2))
        text1 <- c(text1,hits$rsid[i])
        text1 <- c(text1,round(hits$P[i],2))
        
        for(j in 2:3)
        {
          text(x=textpos[j],y=midsplits[i],labels = text1[j], cex = annotext, adj = c(0,0.5))
        }
        gene_info <- ginfo[ginfo$rsid==hits$rsid[i],c("name2","lower","upper")]
        dist <- paste(round(as.numeric(gene_info$upper[1])/1e6,2),
                      round(as.numeric(gene_info$lower[1])/1e6,2), sep = " - ")
        #dist <- paste0(round(as.numeric(dist)/1e4),"Kb")
        text(x=textpos[5],y=midsplits[i],
             labels = paste0(paste0(hits$chrom[i],":",dist),
                             "\n(",round(round(as.numeric(gene_info$lower[1])/1e6,2)-round(as.numeric(gene_info$upper[1])/1e6,2),2),"Mb)"),
             cex = annotext, adj = c(0,0.5), font = 1)
        if(nrow(gene_info) == 1 ) #|| length(unique(gene_info$lower)) == 1)
        {
          dist <- paste(round(as.numeric(gene_info$upper[1])/1e6,2),
                        round(as.numeric(gene_info$lower[1])/1e6,2), sep = "-")
          #dist <- paste0(round(as.numeric(dist)/1e4),"Kb")
          gene_name_plot <- unique(gene_info$name2)
          if(is.na(gene_name_plot)) gene_name_plot <- "No nearby genes"
          text(x=textpos[4],y=midsplits[i],labels = paste(gene_name_plot, collapse = " "),
               cex = annotext, adj = c(0,0.5), font = 3, xpd = T)
          #text(x=textpos[5],y=midsplits[i],labels = dist, cex = 0.8, adj = c(0,0.5), font = 1)
        } else
        {
          num_genes <- nrow(gene_info)
          gene_names <- paste(gene_info$name2,collapse =" ")
          text(x=textpos[5],y=midsplits[i],labels = gene_info$distance[1],
               cex = annotext, adj = c(0,0.5), font = 1, xpd = T)
          ngene_names <- num_gene_rows[(i+1)]
          if(ngene_names %in% c(1,2))
          {
            diff_midsplits <- midsplits
            if(length(midsplits)> 1) diff_midsplits <- diff(midsplits)
            
            newsplits <- seq(midsplits[i]-(diff_midsplits[1]*0.175),midsplits[i]+(diff_midsplits[1]*0.175),length=2)
            spaces <- grep(" ",strsplit(gene_names,"")[[1]])
            nearest <- spaces[which.min(abs((35)-spaces))]
            gene_names <- c(substr(gene_names,1,nearest-1),
                            substr(gene_names,nearest+1,nchar(gene_names)))
            if(nearest<=35)
            {
              text(x=textpos[4],y=midsplits[i],labels = paste(gene_names,collapse=" "),
                   cex = annotext, adj = c(0,0.5), font = 3, xpd = T) 
            } else
            {
              for(j in 1:2)
              {
                text(x=textpos[4],y=newsplits[j],labels = gene_names[j],
                     cex = annotext, adj = c(0,0.5), font = 3, xpd = T)
              }
            }
          }
          if(ngene_names == 3)
          {
            diff_midsplits <- midsplits
            if(length(midsplits)> 1) diff_midsplits <- diff(midsplits)
            newsplits <- seq(midsplits[i]-(diff_midsplits[1]*0.35),midsplits[i]+(diff_midsplits[1]*0.35),length=3)
            spaces <- grep(" ",strsplit(gene_names,"")[[1]])
            nearest1 <- spaces[which.min(abs((35)-spaces))]
            nearest2 <- spaces[which.min(abs((70)-spaces))]
            
            gene_names <- c(substr(gene_names,1,nearest1-1),
                            substr(gene_names,nearest1+1,nearest2-1),
                            substr(gene_names,nearest2+1,nchar(gene_names)))
            for(j in 1:3)
            {
              text(x=textpos[4],y=newsplits[j],labels = gene_names[j],
                   cex = annotext, adj = c(0,0.5), font = 3, xpd = T)
            }
          }
          if(ngene_names == 4)
          {
            diff_midsplits <- midsplits
            if(length(midsplits)> 1) diff_midsplits <- diff(midsplits)
            newsplits <- seq(midsplits[i]-(diff_midsplits[1]*0.35),midsplits[i]+(diff_midsplits[1]*0.35),length=4)
            spaces <- grep(" ",strsplit(gene_names,"")[[1]])
            nearest1 <- spaces[which.min(abs((35)-spaces))]
            nearest2 <- spaces[which.min(abs((70)-spaces))]
            nearest3 <- spaces[which.min(abs((105)-spaces))]
            
            gene_names <- c(substr(gene_names,1,nearest1-1),
                            substr(gene_names,nearest1+1,nearest2-1),
                            substr(gene_names,nearest2+1,nearest3-1),
                            substr(gene_names,nearest3+1,nchar(gene_names)))
            for(j in 1:4)
            {
              text(x=textpos[4],y=newsplits[j],labels = gene_names[j],
                   cex = annotext, adj = c(0,0.5), font = 3, xpd = T)
            }
          }
          if(ngene_names == 5)
          {
            diff_midsplits <- midsplits
            if(length(midsplits)> 1) diff_midsplits <- diff(midsplits)
            newsplits <- seq(midsplits[i]-(diff_midsplits[1]*0.5),midsplits[i]+(diff_midsplits[1]*0.5),length=5)
            spaces <- grep(" ",strsplit(gene_names,"")[[1]])
            nearest1 <- spaces[which.min(abs((35)-spaces))]
            nearest2 <- spaces[which.min(abs((70)-spaces))]
            nearest3 <- spaces[which.min(abs((105)-spaces))]
            nearest4 <- spaces[which.min(abs((140)-spaces))]

            gene_names <- c(substr(gene_names,1,nearest1-1),
                            substr(gene_names,nearest1+1,nearest2-1),
                            substr(gene_names,nearest2+1,nearest3-1),
                            substr(gene_names,nearest3+1,nearest4-1),
                            substr(gene_names,nearest4+1,nchar(gene_names)))
            for(j in 1:5)
            {
              text(x=textpos[4],y=newsplits[j],labels = gene_names[j],
                   cex = annotext, adj = c(0,0.5), font = 3, xpd = T)
            }
          }
        }
      }
        
      for(i in 2:5)
      {
        text(x=textpos[i],y=midsplits[1]-diff(midsplits)[1],
             labels = texttitles[i], cex = 1,adj = c(0,0.5), font = 2, xpd = T, srt = 0)
      }
    }  
    ##############################################################
    ## PLOT BETAS
    par(mar=c(botmar,0,topmar,0))
    xlims <- c(-4,4)
    plot(rep(0,nrow(snps)),pos2,
         ylab="",xlab="",yaxs = "i",xaxs = "i",
         xaxt="n",type="n",yaxt="n",axes = F,
         xlim=xlims,ylim=c(xmin,xmax))
    abline(v=0, xpd = F)
    axis(1)
    
    if(nrow(hits) > 0)
    {
      for(i in 1:nrow(hits))
      {
        beta <- hits$beta[i]
        if(length(midsplits) == 1)
        {
          middiff <- midsplits
           
        } else
        {
          middiff <- diff(midsplits)[1]
        }
        newsplits <- seq(midsplits[i]-(middiff*0.25),midsplits[i]+(middiff*0.25),length=2)
        betacol <- pcolshex[ancreg_list==gsub("\\.","\\-",hits$region[i])]
        abline(h=midsplits[i],lty=2, xpd = F)
        rect(0,newsplits[1],beta,newsplits[2], col = betacol, border = betacol)
      }
      text(x=0,y=midsplits[1]-diff(midsplits)[1],
           labels = expression(Delta~Ancestry), cex = 1,adj = c(0.5,0.5), font = 2, xpd = T,
           srt = 0)
      text(x=-3,y=midsplits[1]-diff(midsplits)[1],
           labels = "-", cex = 1,adj = c(0,0.5), font = 2, xpd = T, srt = 0)
      text(x=3,y=midsplits[1]-diff(midsplits)[1],
           labels = "+", cex = 1,adj = c(0,0.5), font = 2, xpd = T, srt = 0)
    }  
    dev.off()
  }
}


### functions based on HLA-TAPAS manhanttan plot module, tweeking to plot classical alleles next to amino acids

# function itself
make.fancy.locus.plot.bare <- function(chr, title, locus, omnibus, min.pos, max.pos, yrange, hitsnp, r.data, yax, 
                                       p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "#F5BE41", p.color.SNP = "#626262",
                                       topcolor = "#FF0000", cex=cex, topcex=topcex,
                                       hla.wiggle = 3e5, aa.wiggle = 3e5, snp.wiggle = 3e5){
  

  locus <- locus[order(-locus$LOG10P),]
  par(mar=c(0,5.5,3,0))
  plot(0, xlim=c(min.pos, max.pos), ylim=c(0,yrange), xlab="", ylab="",  axes=F)
  title(title, cex.main = cex, font.main =2, adj = 0.01, line = 0)
  axis(2, at=yax, las=1, cex.axis=1.5)
  axis(1, at=seq(28,34,1)*1E6, labels=rep("",7),line=0.2)
  mtext(text=bquote(-log[10]~italic("P")), side=2, at=(yrange/2), line=3.5, cex=1.5)
  
  lines(c(min.pos, max.pos), c(-log10(5E-8),-log10(5E-8)), lty="longdash", lwd=1, col="#333333")
  
  ### `locus` will be divided into 3 groups (2019. 04. 01.)
  # 1. Normal Variants(Intergenic or Intragenic SNPs) : locus.Variants
  # 2. HLA markers : locus.HLA
  # 3. AA markers : locus.AA
  # 4. HLA Intragenic SNPs : locus.SNP
  
  col.label = locus[, 1]
  
  f_AA = grepl("^AA_", col.label)
  f_HLA = grepl("^HLA_", col.label)
  f_SNP = grepl("^SNP(S)?_(.+)_", col.label)
  
  f_HLA_Markers = (f_AA | f_HLA) | f_SNP
  f_Variant = !f_HLA_Markers
  
  locus.AA = locus[f_AA, ]
  locus.HLA = locus[f_HLA, ]
  locus.SNP = locus[f_SNP, ]
  locus.Variant = locus[f_Variant, ]
  
  hitAA <- na.omit(omnibus) %>% filter(LOG10P == min(LOG10P)) %>% select(SNP) %>% as.matrix() %>% as.character() 
  hitHLA <- na.omit(locus.HLA) %>% filter(LOG10P == min(LOG10P)) %>% select(SNP) %>% as.matrix() %>% as.character() 
  hitVariant <- na.omit(locus.Variant) %>% filter(LOG10P == min(LOG10P)) %>% select(SNP) %>% as.matrix() %>% as.character() 
  
  hitHLA <- hitHLA[1]
  
  points(locus.Variant$BP, -locus.Variant$LOG10P, pch=20, cex=cex, lwd=0, bg=p.color.Variant, col=p.color.Variant)
  points(locus.SNP$BP, -locus.SNP$LOG10P, pch=20, cex=cex, lwd=0, bg=p.color.SNP, col=p.color.SNP)
  points(locus.AA$BP, -locus.AA$LOG10P, pch=20, cex=cex, lwd=0, bg=p.color.AA, col=p.color.AA)
  points(omnibus$BP, -omnibus$LOG10P, pch=23, cex=cex, lwd=0, bg=p.color.AA,col=p.color.AA)
  points(locus.HLA$BP, -locus.HLA$LOG10P, pch=20, cex=cex, lwd=0, bg=p.color.HLA, col=p.color.HLA)
  
  if (hitAA != "") {
    hit <- omnibus[omnibus$SNP==hitAA,]
    points(hit$BP, -hit$LOG10P, pch=23, cex=topcex, lwd=1.5, col=p.color.AA)
    labelAA <- paste(unlist(strsplit(hitAA,"_"))[1:3],collapse = "_")
    text(hit$BP + aa.wiggle, -hit$LOG10P,labels=labelAA, col = p.color.AA,cex=cex)
  }
  
  if (hitHLA != "") {
    hit <- locus[locus$SNP==hitHLA,]
    points(hit$BP, -hit$LOG10P, pch=20, cex=topcex, lwd=1.5, col=p.color.HLA)
    if (length(unlist(strsplit(hitHLA,":"))) > 2){
      labelHLA <- paste(unlist(strsplit(hitHLA,":"))[1:2],collapse = ":")
    }
    else
      labelHLA <- hitHLA
      
    text(hit$BP + hla.wiggle, -hit$LOG10P,labels=labelHLA, col = p.color.HLA,cex=cex)
  }
  
  if (hitVariant != "") {
    hit <- locus[locus$SNP==hitVariant,]
    points(hit$BP, -hit$LOG10P, pch=20, cex=topcex, lwd=1.5, col=p.color.Variant)
    text(hit$BP + snp.wiggle, -hit$LOG10P,labels=hitVariant, col = p.color.Variant,cex=cex)
  }
}


make.fancy.locus.plot.bottom <- function(chr, min.pos, max.pos, pathToTheGeneBuild) {
  ## Plot margins
  par(mar=c(0,5.5,0.5,0))
  #par(mar=c(4,5.5,0.5,4))
  ##
  ## yrange of y-axis
  ##
  ## this dedicates proportion of the yaxis to the genes, labels, recomb rate
  yrange = 0
  offset <- 10
  big.yrange <- yrange + offset
  ystart.recomb <- -yrange/10
  gene.y <- -offset * 0.4 + ystart.recomb * 0.4
  plot(-1, -1, xlim=c(min.pos, max.pos), ylim=c(-offset,yrange), xlab="", ylab="", axes=F)
  
  ##
  ## genes in the region (build 36)
  ##
  
  genelist <- read.table(pathToTheGeneBuild, header=T)
  
  # ### (2018. 9. 24.) Modified to work with pure "knownGene.txt" file distributed by UCSC Annotation Database.
  # genelist <- read.table(pathToTheGeneBuild, header=F, sep = '\t')
  # colnames(genelist) = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  
  # # colnames for "knownGene"  
  # c("name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID")
  
  genelist <- subset(genelist, genelist$chrom == paste("chr", chr, sep=""))
  genes.in.locus <- subset(genelist, (genelist$txStart > min.pos & genelist$txStart < max.pos ) | (genelist$txEnd > min.pos & genelist$txEnd < max.pos))
  genes.in.locus <- genes.in.locus[,c("txStart", "txEnd", "name2")]
  names(genes.in.locus) <- c("START", "END", "GENE")
  
  
  ##
  ## axes, titles and legends
  ##
  mtext(paste("Chromosome", chr, "position (Mb)", sep=" "), side=1, line=2.5, cex=1.5)
  axis(1, at=seq(28,34,1)*1E6, labels=seq(28,34,1), cex.axis=1.5)
  
  ##
  ## plot the genes
  ##
  boxheight = big.yrange / 4
  geneboxheight = boxheight * 1.3
  for ( i in 1:nrow(genes.in.locus) ) {
    lines(rep(mean(c(genes.in.locus[i,]$START,genes.in.locus[i,]$STOP)),2), c(gene.y+boxheight/2,gene.y-boxheight/2), lwd=1.5, col="darkgrey")
  }
  j=1
  for ( i in 1:nrow(genes.in.locus) ) {
    genename = genes.in.locus[i,]$GENE
    # if (genename %in% strsplit("HLA-A HLA-C HLA-B TNF HLA-DRA HLA-DRB1 HLA-DQA1 HLA-DQB1 HLA-DPA1 HLA-DPB1", " ")[[1]]) {
    if (genename %in% strsplit("HLA-A HLA-C HLA-B HLA-DRB1 HLA-DQA1 HLA-DQB1 HLA-DPA1 HLA-DPB1", " ")[[1]]) {
      genecols = strsplit("#00ab34 #c800a0 #7100c8 #c80000 #567200 #ecbd00 #FF0000 #0033F2"," ")[[1]]
      # genecols = strsplit("#BF00EF black #00C418 black black #FF0000 black black black #0033F2"," ")[[1]]
      
      ## genecols = rep("black", 8)
      genetextcols = rep("black", 20)
      
      # genetips = c("upleft", "upleft", "upright", "downleft", "downleft", "upleft", "upright", "downright", "upright", "downright")
      genetips = c("upleft", "upleft", "upright", "upleft", "upright", "downright", "upright", "downright")
      tiplen.y = boxheight / 2
      tiplen.x = (max.pos-min.pos)/75
      pickcol = genecols[j]
      picktextcol = genetextcols[j]
      gene.center = mean(c(genes.in.locus[i,]$START,genes.in.locus[i,]$STOP))
      gene.head = gene.y+geneboxheight/2
      gene.foot = gene.y-geneboxheight/2
      lines(rep(gene.center,2), c(gene.head, gene.foot), lwd=3, col=pickcol)
      k = 1.7
      textcex = 2
      textoffset = 0
      if (genename == "HLA-DQA1" || genename == "HLA-DQB1") {
        textoffset = -1 ## To avoid overlap.
      }
      if (length(grep("up",genetips[j])) > 0) {
        if (length(grep("left",genetips[j])) > 0) {
          lines(c(gene.center, gene.center - tiplen.x), c(gene.head, gene.head + tiplen.y), lwd=0.4, col="black")
          text(x=gene.center - tiplen.x, gene.head + k*tiplen.y, label=genename,
               col=picktextcol, cex=textcex, font=3, pos=2, xpd=T, offset=textoffset)
        } else {
          lines(c(gene.center, gene.center + tiplen.x), c(gene.head, gene.head + tiplen.y), lwd=0.4, col="black")
          text(x=gene.center + tiplen.x, gene.head + k*tiplen.y, label=genename,
               col=picktextcol, cex=textcex, font=3, pos=4, xpd=T, offset=textoffset)
        }
      } else {
        if (length(grep("left",genetips[j])) > 0) {
          lines(c(gene.center, gene.center - tiplen.x), c(gene.foot, gene.foot - tiplen.y), lwd=0.4, col="black")
          text(x=gene.center - tiplen.x, gene.foot - k*tiplen.y*1.2, label=genename,
               col=picktextcol, cex=textcex, font=3, pos=2, xpd=T, offset=textoffset)
        } else {
          lines(c(gene.center, gene.center + tiplen.x), c(gene.foot, gene.foot - tiplen.y), lwd=0.4, col="black")
          text(x=gene.center + tiplen.x, gene.foot - k*tiplen.y*1.2, label=genename,
               col=picktextcol, cex=textcex, font=3, pos=4, xpd=T, offset=textoffset)
        }
      }
      print(c(toString(genes.in.locus[i,]$GENE),
              toString(genes.in.locus[i,]$START),
              toString(genes.in.locus[i,]$END),
              toString(genes.in.locus[i,]$END-genes.in.locus[i,]$START),
              toString(pickcol)))
      j=j+1
    }
  }
}

TRIM_LOGISTIC.ASSOC2 = function(p_assoc.logistic_,bim=bim){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  HLA.str <- df_assoc.logistic[grepl("^HLA_",df_assoc.logistic$ID),]$ID
  HLA.ID <- sapply(HLA.str,function(x) {str <- unlist(strsplit(x,"\\.")); allele<-paste(str[-1],collapse = ":"); paste(str[1],"*",allele,sep="")})
  df_assoc.logistic$ID[grepl("^HLA_",df_assoc.logistic$ID)] <- as.character(HLA.ID)
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  out <- cbind(df_assoc.logistic[,c("ID", "BP", "P")], LOG10P)
  names(out) <- c("SNP","BP","P","LOG10P")
  return(out)
}

TRIM_LOGISTIC.ASSOC3 = function(p_assoc.logistic_,bim=bim){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  out <- cbind(df_assoc.logistic[,c("ID", "BP", "P")], LOG10P)
  names(out) <- c("SNP","BP","P","LOG10P")
  return(out)
}

TRIM_LOGISTIC.ASSOC = function(p_assoc.logistic_){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  out <- cbind(df_assoc.logistic[,c("SNP", "BP", "P")], LOG10P)
  return(out)
}




TRIM_OMNI.ASSOC <- function(p_omnibus){
  df_assoc.omnibus = read.table(p_omnibus, header = T)
  df_omnibus = df_assoc.omnibus %>% filter(AA_POS > 0 )
  LOG10P = log10(df_omnibus$PVALUE)
  out <- cbind(df_omnibus[,c("AA_ID","POS","PVALUE")],LOG10P)
  names(out) <- c("SNP","BP","P","LOG10P")
  return(out)
}

TRIM.META = function(p_assoc.logistic_,bim=bim){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  HLA.str <- df_assoc.logistic[grepl("^HLA_",df_assoc.logistic$ID),]$ID
  HLA.ID <- sapply(HLA.str,function(x) {str <- unlist(strsplit(x,"\\.")); allele<-paste(str[-1],collapse = ":"); paste(str[1],"*",allele,sep="")})
  df_assoc.logistic$ID[grepl("^HLA_",df_assoc.logistic$ID)] <- as.character(HLA.ID)
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  df_assoc.logistic$REF <- bim[match(df_assoc.logistic$ID,bim$V2),5]
  df_assoc.logistic$ALT <- bim[match(df_assoc.logistic$ID,bim$V2),6]
  df_assoc.logistic[grepl("^HLA_",df_assoc.logistic$ID),]$REF <- "Present"
  df_assoc.logistic[grepl("^HLA_",df_assoc.logistic$ID),]$REF <- "Absent"
  out <- df_assoc.logistic[,c("ID", "BP","REF","ALT","BETA","SE", "P")]
  return(out)
}

TRIM.META2 = function(p_assoc.logistic_,bim=bim,freq=frq){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  df_assoc.logistic$REF <- bim[match(df_assoc.logistic$ID,bim$V2),5]
  df_assoc.logistic$ALT <- bim[match(df_assoc.logistic$ID,bim$V2),6]
  df_assoc.logistic$MAF <- frq[match(df_assoc.logistic$ID,frq$SNP),]$MAF
  out <- df_assoc.logistic[,c("ID", "BP","REF","ALT","MAF","BETA","SE", "P")]
  return(out)
}


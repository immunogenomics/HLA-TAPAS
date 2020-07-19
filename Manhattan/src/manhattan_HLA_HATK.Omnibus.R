########## <Argument Pasrsing> ##########
args = commandArgs(trailingOnly = TRUE)


p_assoc.logsitic_ = strsplit(args[1], ',')[[1]]
out_ = args[2]

# pointcol_ = args[]
pointsize_ = as.numeric(args[3])
topcol_ = args[4]
min.pos_ = as.numeric(args[5])
max.pos_ = as.numeric(args[6])
topsignal_ = strsplit(args[7], ',')[[1]] # label of top signal marker
max.yaxis_ = as.numeric(strsplit(args[8], ',')[[1]])
yaxis_unit_ = as.numeric(args[9])
p_src_ = args[10]

cat("<Given arguments from Omnibus Manhattan.>\n")
print(args)



make.fancy.locus.plot.bare <- function(chr, title, locus, min.pos, max.pos, yrange, hitsnp, r.data, yax, 
                                       p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "#F5BE41", p.color.SNP = "#626262",
                                       topcolor = "#FF0000", arrowsnp="", hitsnpB="", hitsnpC="") {
  
  
  par(mar=c(0,5.5,2,4))
  plot(0, type="n", xlim=c(min(0,min.pos), max.pos), ylim=c(0,yrange), xlab="", ylab="", main=title, axes=F)
  axis(1, at=if(min.pos < 0) c(min.pos, 0, max.pos) else c(min.pos, max.pos), labels= if(min.pos < 0) c(min.pos, 0, max.pos) else c(min.pos, max.pos), line=0.2)
  axis(2, at=yax, las=1, cex.axis=1.5)
  mtext(text=bquote(-log[10]~italic("P")), side=2, at=(yrange/2), line=3.5, cex=1.7)
  lines(c(min.pos, max.pos), c(-log10(5E-8),-log10(5E-8)), lty="longdash", lwd=1, col="#333333")
  points(locus$REL_POS, -locus$LOG10P, pch=23, cex=2.3, lwd=0, bg=p.color.AA, col=p.color.AA)
  
  
  # (Top signal)
  # if (hitsnp != "") {
  #     hit <- locus[locus$SNP==hitsnp,]
  #     points(hit$BP, -hit$LOG10P, pch=23, cex=2.4, lwd=1.5, bg=topcolor)
  #     text(x = hit$BP, y = -hit$LOG10P, labels = hitsnp, font = 2, family="sans", adj = c(NA, -2.0), cex = 1.0)

  # }

  if (hitsnpB != "") {
     hit <- locus[locus$SNP==hitsnpB,]
     points(hit$BP, -hit$LOG10P, pch=23, cex=2.4, lwd=1.5, bg=topcolor)
  }

  if (hitsnpC != "") {
     hit <- locus[locus$SNP==hitsnpC,]
     points(hit$BP, -hit$LOG10P, pch=23, cex=2.4, lwd=1.5, bg=topcolor)
  }


  if (arrowsnp != "") {
      arro <- locus[locus$RSID==arrowsnp,]
      height <- max(-locus[locus$BP < arro$BP+50000 & locus$BP > arro$BP-50000,]$LOG10P)
      arrows(x0=arro$BP, y0=height+yrange/5, y1=height+yrange/20, length=0.25, code=2, col="red", lwd=5)
  }
}

########## <Loading Necessary Scripts(Modules)> ##########


#### Main

pdf(paste0(out_, ".pdf"), width=16, height=9, pointsize=pointsize_) # `pointsize` is now manipulated by the argument 'pointsize_'. (2019. 04. 09.)

NumberofLogistic = length(p_assoc.logsitic_)
layout(matrix(1:(NumberofLogistic+1), (NumberofLogistic+1),1, byrow=TRUE), heights=c(rep(7, NumberofLogistic), 2))
# (2019. 04. 01.) One manhattan graph : One fancy buttom = 7 : 2

# (2018. 8. 10.) No more `r.data` for colorRamp().

for (i in 1:NumberofLogistic) {

  cat('\n(', i, ')\n')
  print(p_assoc.logsitic_[i])
  print(topsignal_[i])
  print(max.yaxis_[i])

  locus <- read.table(p_assoc.logsitic_, header = T)
  locus = locus[, c('Variant', 'REL_POS', 'log10_P')]
  colnames(locus) = c('Variant', 'REL_POS', 'LOG10P')
  # manhattan plot
  make.fancy.locus.plot.bare("6", (if(i == 1) "" else ""), locus, min.pos_, max.pos_, max.yaxis_[i], topsignal_[i], NULL, seq(0,max.yaxis_[i],by=yaxis_unit_), topcolor=topcol_) #
  
}

dev.off()

cat("<Manhattan Plotting done.>\n")
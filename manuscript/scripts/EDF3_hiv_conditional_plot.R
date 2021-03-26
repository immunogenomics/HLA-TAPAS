library(dplyr)
source("MHC_locus_plot.R")
#setwd("~/GitHub/raychaudhuri-lab/HLA-TAPAS//manuscript/scripts/")
bim <- read.table("../data/all.bim",h=F,stringsAsFactors = F)
frq <- read.table("../data/all.frq",h=T,stringsAsFactors = F)
min.pos <- 28e6
max.pos <- 34e6
cex <- 2
topcex <- 3

out <- "../figs/EDF3_HIV_conditional"
pdf(paste0(out, ".pdf"),width=18,height=20,res=600) 
#png(paste0(out, ".png"), width=18, height=20,res=200,unit="in")

NumberofLogistic = 5
layout(matrix(1:(NumberofLogistic+1), (NumberofLogistic+1),1, byrow=TRUE), heights=c(7,6,5,4,4, 2))

locus <- TRIM_LOGISTIC.ASSOC3("../data/all.1st.assoc.linear",bim=bim)
omnibus <- TRIM_OMNI.ASSOC("../data/all.1st.txt")
make.fancy.locus.plot.bare(chr = "6",locus=locus, omnibus = omnibus,
                           title = "(a) no condition",min.pos = min.pos,max.pos = max.pos, yrange=200,
                           "",yax = NULL, cex = cex, topcex = topcex, hla.wiggle = -5e5, snp.wiggle = 7e5,aa.wiggle = 5e5,
                           p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")


locus <- TRIM_LOGISTIC.ASSOC2("../data/all.2nd.assoc.linear",bim=bim)
omnibus <- TRIM_OMNI.ASSOC("../data/all.2nd.txt")
make.fancy.locus.plot.bare(chr = "6",locus=locus %>% filter(locus$SNP %in% frq[frq$MAF>0.005,]$SNP), omnibus = omnibus,
                           title = "(b) control for position 97 in HLA-B",min.pos = min.pos,max.pos = max.pos, yrange=45,
                           "",yax = NULL, cex = cex, topcex = topcex,snp.wiggle = -5e5,aa.wiggle = 5e5,hla.wiggle = 6e5,
                           p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")

locus <- TRIM_LOGISTIC.ASSOC2("../data/all.3rd.assoc.linear",bim=bim)
omnibus <- TRIM_OMNI.ASSOC("../data/all.3rd.txt")
make.fancy.locus.plot.bare(chr = "6",locus=locus %>% filter(locus$SNP %in% frq[frq$MAF>0.005,]$SNP), omnibus = omnibus,
                           title = "(c) control for position 97 and 67 in HLA-B",min.pos = min.pos,max.pos = max.pos, yrange=35,
                           "",yax = NULL, cex = cex, topcex = topcex, aa.wiggle = -6e5,hla.wiggle = 5e5,snp.wiggle = 6e5,
                           p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")


locus <- TRIM_LOGISTIC.ASSOC2("../data/all.4th.assoc.linear",bim=bim)
omnibus <- TRIM_OMNI.ASSOC("../data/all.4th.txt")
make.fancy.locus.plot.bare(chr = "6",locus=locus %>% filter(locus$SNP %in% frq[frq$MAF>0.005,]$SNP), omnibus = omnibus,
                           title = "(d) control for position 97, 67 and 156 in HLA-B",min.pos = min.pos,max.pos = max.pos, yrange=20,
                           "",yax = NULL, cex = cex, topcex = topcex,aa.wiggle = -5e5,hla.wiggle = 5e5,snp.wiggle = 6e5,
                           p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")

locus <- TRIM_LOGISTIC.ASSOC2("../data/all.allB.assoc.linear",bim=bim)
omnibus <- TRIM_OMNI.ASSOC("../data/all.allB.txt.gz")
make.fancy.locus.plot.bare(chr = "6",locus=locus, omnibus = omnibus,
                           title = "(e) control for all amino acid positions in HLA-B",min.pos = min.pos,max.pos = max.pos, yrange=20,
                           "",yax = NULL, cex = cex, topcex = topcex,aa.wiggle = -5e5,snp.wiggle = 6e5,hla.wiggle = 5e5,
                           p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")


make.fancy.locus.plot.bottom("6", min.pos, max.pos, pathToTheGeneBuild="~/GitHub/yang/HLA-TAPAS/Manhattan/data/known_genes/known_genes_chr6.hg19.txt")
dev.off()

# pdf(paste0(out, ".pdf"),paper="a4") 
# # png(paste0(out, ".png"), width=24, height=16,res=200,unit="in") 
#  layout(matrix(1:(NumberofLogistic+1), (NumberofLogistic+1),1, byrow=TRUE), heights=c(6, 2))
# make.fancy.locus.plot.bare(chr = "6",locus=locus, omnibus = omnibus,
#  title = "",min.pos = min.pos,max.pos = max.pos, yrange=200,
#                 "",yax = NULL, cex = cex, topcex = topcex, hla.wiggle = 4e5, snp.wiggle = -5e5,
#  p.color.Variant = "#999999", p.color.AA = "#CF3721", p.color.HLA = "darkblue", p.color.SNP = "#999999")
#  
# make.fancy.locus.plot.bottom("6", min.pos, max.pos, pathToTheGeneBuild="~/GitHub/yang/HLA-TAPAS/Manhattan/data/known_genes/known_genes_chr6.hg19.txt")
# dev.off()

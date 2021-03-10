# code for regenerating Figure 4 in Luo et al. 2020

library(ggrepel)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)

source("amino_acid_manhattan_plot.R")

# set global plotting parameters
theme_set(theme_classic() +
            theme(text = element_text(size=20, family = "sans")))

pop_colors <- c("#938dd2","#E69F00" ,"#56B4E9","#D55E00","#009E73")
pop_labels <- c("Admixed African","East Asian","European","Latino", "South Asian")

#Calc color palette (discrete)
gene_colors <- c("A"="#004586", "B" = "#ff420e", "C"="#ffd320", 
                 "DQA1" = "#579d1c", "DQB1" = "#7e0021", "DRB1" = "#83caff","DPA1" = "#314004", "DPB1" = "#aecf00"  )




bim <- read.table("../data/all.bim",h=F,stringsAsFactors = F)
df <- read.table("../data/all.1st.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p1 <- aa.manhattan(assoc, showxlab=FALSE,hlight=c("B-97"),ymax=210)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC("../data/all.1st.assoc.linear",bim=bim)
p2 <- hla.manhattan(df = hla, showxlab=FALSE,ymax=210,nudge_value = 10)

# 2nd round
df <- read.table("../data/all.2nd.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p3 <- aa.manhattan(assoc, showxlab=FALSE,hlight=c("B-67"),ymax=50)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC2("../data/all.2nd.assoc.linear",bim=bim)
p4 <- hla.manhattan(df = hla, showxlab=FALSE,ymax=50,hlight=c("B*81:01:01G"))

# 3rd round
df <- read.table("../data/all.3rd.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p5 <- aa.manhattan(assoc, showxlab=FALSE,hlight=c("B-156"),ymax=35)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC2("../data/all.3rd.assoc.linear",bim=bim)
p6 <- hla.manhattan(df = hla, showxlab=FALSE,ymax=35,hlight=c("B*81:01:01G"))

# 4th round
df <- read.table("../data/all.4th.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p7 <- aa.manhattan(assoc, showxlab=TRUE,hlight=c("A-77"),ymax=10)

hla <- TRIM_LOGISTIC.HLA.ASSOC2("../data/all.allB.assoc.linear",bim=bim)
p8 <- hla.manhattan(df = hla, showxlab=TRUE,ymax=10, nudge_value = 1)


png("../figs/Fig4_conditional_assoc_withHLA.png", height = 13, width = 12, units= "in", res=200)
y.grob <- textGrob("-log10(P)", 
                   gp=gpar(fontface="bold", col="black", fontsize=20), rot=90, vjust = .3)

p <- plot_grid(p2, p1,p4,p3, p6,p5,p8,p7, align = "v", nrow = 4, rel_widths = c(1/4, 3/4), 
               rel_heights = c(1.8/4, 1.4/4, 1/4,.8/4),
               labels = c("(a)","(b)", "(c)","(d)", "(e)","(f)", "(g)","(h)"),  hjust = .3, vjust = 1)

grid.arrange(arrangeGrob(p, left = y.grob))

dev.off()



df <- read.table("../data/all.1st.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p1 <- aa.manhattan(assoc, showxlab=FALSE,hlight=c("B-97"),ymax=210)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC("../data/all.1st.assoc.linear",bim=bim)
p2 <- hla.manhattan(df = hla, showxlab=FALSE,ymax=210,nudge_value = 10)

# 2nd round
df <- read.table("../data/all.2nd.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p3 <- aa.manhattan(assoc, showxlab=FALSE,hlight=c("B-67"),ymax=50)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC2("../data/all.2nd.assoc.linear",bim=bim)
p4 <- hla.manhattan(df = hla, showxlab=FALSE,ymax=50,hlight=c("B*81:01:01G"))

# 3rd round
df <- read.table("../data/all.3rd.txt",h=T,stringsAsFactors = F)
assoc <- df[df$AA_POS>0,]
assoc$SNP <- paste(assoc$GENE,assoc$AA_POS,sep="-")
p5 <- aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-156"),ymax=35)
#ggsave(aa.manhattan(assoc, showxlab=TRUE,hlight=c("B-97")),filename = "../figs/aa_manhattan_97.png",height=7,width=12)

hla <- TRIM_LOGISTIC.HLA.ASSOC2("../data/all.3rd.assoc.linear",bim=bim)
p6 <- hla.manhattan(df = hla, showxlab=TRUE,ymax=35,hlight=c("B*81:01:01G"))

png("~/presentations/jobs/hiv_aa_cond2.png", height = 7, width = 12, units= "in", res=200)
y.grob <- textGrob("-log10(P)", 
                   gp=gpar(fontface="bold", col="black", fontsize=20), rot=90, vjust = .3)



p <- plot_grid(p2, p1,p4,p3,p6,p5, align = "v", nrow = 3, rel_widths = c(1/4, 3/4), 
               rel_heights = c(1/2, 1/4,1/4),
               labels = c("(a)","(b)", "(c)","(d)", "(e)","(f)"),  hjust = .3, vjust = 1)

grid.arrange(arrangeGrob(p, left = y.grob))

dev.off()

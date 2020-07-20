# scale by amino acid
aa.manhattan <- function(df = assoc, showxlab = TRUE,  showylab = FALSE, hlight=NULL,
                         nudge_value = -.6, ymax=NULL){
  df.assoc <- NULL
  genes <- c("A","C","B","DRB1","DQA1","DQB1","DPA1","DPB1")
  for (gene in genes){
    assoc <- df %>% filter(grepl("AA_",AA_ID)) %>% na.omit() %>% 
      select(GENE,AA_POS,POS,PVALUE, SNP) %>%  group_by(GENE) %>% arrange(AA_POS)
    df1 <- assoc
    df.tmp <- df1[df1$GENE==gene, ]
    df.tmp$BP <- (df.tmp$AA_POS - df.tmp$AA_POS[1])*3
    df.assoc <- rbind(df.assoc, df.tmp)
  }
  
  # proportional to gene length
  df.assoc$GENE <- factor(df.assoc$GENE, levels=genes)
  df.assoc$genenum <- as.numeric(factor(df.assoc$GENE))
  
  df.tmp <- df.assoc %>%   
    # compute gene length
    group_by(GENE,POS) %>% mutate(rowID = 1:n()) 
  
  
  df.assoc <- df.tmp %>% group_by(GENE) %>%
    summarise(gene_len=max(BP) ) %>% 
    # calcuate cumulative positoin of each gene
    mutate(tot = cumsum(gene_len) - gene_len ) %>%
    
    # add info
    left_join(df.tmp,., by=c("GENE"="GENE")) %>%
    
    # Add a cumulate position of each amino acid
    mutate(BPcum = (tot + rowID)) %>%
    #arrange(GENE) %>% mutate(BPcum = BP/gene_len + genenum * 1.5) %>%
    arrange(GENE,rowID) %>% mutate(BPcum = (BP/gene_len + genenum *1.5)) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% hlight, "yes", "no")) 
  
  
  axisdf <- df.assoc %>% group_by(GENE) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) 
  
  p <- ggplot(df.assoc, aes(x = BPcum, y = -log10(PVALUE), size = -log10(PVALUE), fill=GENE,
                            alpha = -log10(PVALUE))) + 
    geom_point( shape = 23)
  
  if (is.null(ymax)) ymax <- max(-log10(df.assoc$PVALUE))*1.2
  p <- p + scale_fill_manual( values = gene_colors) + 
    scale_x_continuous( label=genes, breaks= axisdf$center ) +
    labs(x = "Amino acid positions") + 
    scale_y_continuous(limits=c(0,ymax), breaks = scales::pretty_breaks(n=3))
  
  # Add highlighted points
  # p <- p + geom_point(data=subset(df.assoc, is_highlight_aa=="yes"), fill="#E64B35FF",size =5, pch=23) +
  # geom_point(data=subset(df.assoc, is_highlight_eur=="yes"), fill="#00A087FF", size = 5, pch=23) 
  p <- p + 
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=df.assoc[df.assoc$is_annotate=="yes",], 
                    aes(label=as.factor(SNP)), size=6,
                    nudge_x = nudge_value, nudge_y = 2, segment.size = 0.2, show.legend = FALSE) 
  
  # add genome-wide threshold
  p <- p + geom_hline(yintercept = -log10(5e-8),linetype="dashed", color="black")
  # whether needs the xlabs
  if(showxlab == TRUE){
    p <- p + theme(legend.position = "none",axis.title.x = element_text(size=14,face="bold",family="sans"),
                   axis.text.x = element_text(angle=90, vjust=.5,hjust=1, size =12, family="sans")) 
    
  }else{
    p<- p + theme(legend.position = "none",axis.title.x = element_blank(), axis.text.x = element_blank())
    
  }
  if(showylab == TRUE){
    p <- p + theme(legend.position = "none") + ylab("-log10(P)")
  }else{
    p <- p + theme(legend.position = "none",axis.title.y = element_blank()) 
    
  }
}


hla.manhattan <- function(df = assoc, showxlab = TRUE,  showylab = FALSE,
                          nudge_value = 4, ymax=NULL,hlight.text=NULL){
  df.assoc <- NULL
  
  genes <- c("A","C","B","DRB1","DQA1","DQB1","DPA1","DPB1")
  for (gene in genes){
    df.tmp <- df[df$GENE==gene, ]
    #df.tmp$BP <- (df.tmp$AA_POS - df.tmp$AA_POS[1])*3
    df.tmp$BP <- abs(df.tmp$POS - df.tmp$POS[1]) # scale by basepair position
    df.assoc <- rbind(df.assoc, df.tmp)
  }
  
  df.tmp <- df.assoc %>%   
    # compute gene length
    group_by(GENE,POS) %>% mutate(rowID = 1:n()) 
  
  # proportional to gene length
  df.assoc$GENE <- factor(df.assoc$GENE, levels=genes)
  df.assoc$genenum <- as.numeric(factor(df.assoc$GENE))
  #df.assoc$TYPE <- ifelse(grepl("^HLA_",df.assoc$SNP),"HLA","AA")
  
  df.tmp <- df.assoc
  
  # hlighting top variants
  #hlight <-df.assoc %>% group_by(TYPE) %>% slice(which.min(PVALUE)) %>% select(SNP)
  
  hlight <- df.assoc[which(df.assoc$PVALUE==min(df.assoc$PVALUE)),]$SNP 
  hlight <- hlight[length(hlight)]
  
  
  df.assoc <- df.tmp %>% group_by(GENE) %>%
    #summarise(gene_len = max(abs(BP))) %>% #summarise by BP
    summarise(gene_len=max(BP) ) %>% #summaryise by length
    #arrange(genes) %>%
    # calcuate cumulative positoin of each gene
    mutate(tot = cumsum(gene_len) - gene_len ) %>%
    
    # add info
    left_join(df.tmp,., by=c("GENE"="GENE")) %>%
    
    # Add a cumulate position of each amino acid
    #mutate(BPcum = (tot + rowID)) %>%
    arrange(GENE) %>% mutate(BPcum = BP/gene_len + genenum * 1.5) %>%
    #arrange(GENE,rowID) %>% mutate(BPcum = (tot + BP)) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% hlight, "yes", "no")) 
  
  axisdf <- df.assoc %>% group_by(GENE) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) 
  
  p <- ggplot(df.assoc, aes(x = BPcum, y = -log10(PVALUE), size = -log10(PVALUE), color=GENE,
                            alpha = -log10(PVALUE)), shape=17) + 
    geom_point()
  
  if (is.null(ymax)) ymax <- max(-log10(df.assoc$PVALUE))*1.2
  p <- p + scale_color_manual( values = gene_colors) + 
    scale_x_continuous( label=genes, breaks= axisdf$center ) +
    labs(x = "Classical allele positions") + 
    scale_y_continuous(limits=c(0,ymax), breaks = scales::pretty_breaks(n=3))
  
  # Add highlighted points
  p <- p + 
    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data=df.assoc[df.assoc$is_annotate=="yes",], 
                    aes(label=ifelse(is.null(hlight.text),unlist(strsplit(as.character(SNP),"_"))[2],
                                     hlight.text)), size=6,
                    nudge_y = nudge_value, segment.size = 0.2, show.legend = FALSE) 
  
  # add genome-wide threshold
  p <- p + geom_hline(yintercept = -log10(5e-8),linetype="dashed", color="black")
  # whether needs the xlabs
  if(showxlab == TRUE){
    p <- p + theme(legend.position = "none",axis.title.x = element_text(size=14,face="bold",family="sans"),
                   axis.text.x = element_text(angle=90, vjust=.5,hjust=1, size =12, family="sans")) 
  }else{
    p<- p + theme(legend.position = "none",axis.title.x = element_blank(),
                  axis.text.x = element_blank())
    
  }
  if(showylab == TRUE){
    p <- p + theme(legend.position = "none") + ylab("-log10(P)")
  }else{
    p <- p + theme(legend.position = "none",axis.title.y = element_blank()) 
    
  }
}

TRIM_LOGISTIC.HLA.ASSOC = function(p_assoc.logistic_,bim=bim){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  out <- df_assoc.logistic[,c("ID", "BP", "P")]
  
  col.label = out[, 1]
  f_HLA = grepl("^HLA_", col.label)
  out.HLA = out[f_HLA, ]
  
  names(out.HLA) <- c("SNP","POS","PVALUE")
  out.HLA$GENE <- sapply(out.HLA$SNP,function(x) unlist(strsplit(unlist(strsplit(x,"\\*"))[1],"_"))[2])
  return(na.omit(out.HLA))
}

TRIM_LOGISTIC.HLA.ASSOC2 = function(p_assoc.logistic_,bim=bim){
  df_assoc.logistic = read.table(p_assoc.logistic_, header = T)
  LOG10P = log10(df_assoc.logistic$P)
  df_assoc.logistic$ID <- sapply(df_assoc.logistic$SNP, function(x){ str <-strsplit(as.character(x), "_")[[1]]; paste(str[-length(str)],collapse="_")})
  HLA.str <- df_assoc.logistic[grepl("^HLA_",df_assoc.logistic$ID),]$ID
  HLA.ID <- sapply(HLA.str,function(x) {str <- unlist(strsplit(x,"\\.")); allele<-paste(str[-1],collapse = ":"); paste(str[1],"*",allele,sep="")})
  df_assoc.logistic$ID[grepl("^HLA_",df_assoc.logistic$ID)] <- as.character(HLA.ID)
  df_assoc.logistic$BP <- bim[match(df_assoc.logistic$ID,bim$V2),4]
  out <- df_assoc.logistic[,c("ID", "BP", "P")]
  
  col.label = out[, 1]
  f_HLA = grepl("^HLA_", col.label)
  out.HLA = out[f_HLA, ]
  
  names(out.HLA) <- c("SNP","POS","PVALUE")
  out.HLA$GENE <- sapply(out.HLA$SNP,function(x) unlist(strsplit(unlist(strsplit(x,"\\*"))[1],"_"))[2])
  return(na.omit(out.HLA))
  
}

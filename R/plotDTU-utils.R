getTotalCount <- function(countData, tx2gene){
  # get tx2gene in better format
  geneForEachTx <- tx2gene$gene_id[match(rownames(countData),tx2gene$isoform_id)]
  geneForEachTx <- as.character(geneForEachTx)
  stopifnot(class(geneForEachTx) %in% c("character", "factor"))
  stopifnot(length(geneForEachTx) == nrow(countData))
  
  # adapted from DEXSeq source code
  forCycle <- split(1:nrow(countData), as.character(geneForEachTx))
  all <- lapply(forCycle, function(i) {
    sct <- countData[i, , drop = FALSE]
    rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[, , drop = FALSE])))
    rownames(rs) <- rownames(sct)
    rs
  })
  
  totalCount <- do.call(rbind, all)
  totalCount <- totalCount[rownames(countData), ]
  return(totalCount)
}

label_facet <- function(txID, padj){
  lev <- levels(as.factor(txID))
  lab <- paste0("empFDR = ", padj)
  names(lab) <- lev
  return(lab)
}

visualize_DTU <- function(object, topTable, contrast, coefficients, groups, summaryStat, tx2gene, transcripts) {
  
  plotList <- list()
  countData <- assay(object)
  
  # To calculate usages, we also need the info of transcripts that we don't want to visualize, i.e. all transcripts of all corresponding genes
  genes <- unique(tx2gene[match(transcripts,tx2gene$isoform_id), "gene_id"])
  
  # Required to keep transcripts in the requested order (and not necessarily the order by which they are in tx2gene or alphabetical)
  transcripts_all <- sapply(genes, function(i){
    as.character(tx2gene[which(tx2gene$gene_id == i), "isoform_id"])
  },simplify = T)
  transcripts_all <- unname(unlist(transcripts_all))
  
  # assign cells to a certain group (violin)
  names(groups) <- c(1:length(groups))
  cell_to_group <- unlist(groups)
  names(cell_to_group) <- rep(paste0("violin", names(groups)),times=lengths(groups))
  
  countData <- countData[transcripts_all,unlist(groups)]
  tx2gene <- tx2gene[tx2gene$isoform_id %in% transcripts_all,]
  
  # get total counts (gene-level count) and usages (relative proportion of expression)
  totalCount <- getTotalCount(countData, tx2gene)
  usage <- countData/totalCount
  
  # subset to only the transcripts we want to plot
  usage <- usage[transcripts,,drop=FALSE]
  totalCount <- totalCount[transcripts,,drop=FALSE]
  
  # plot each of the requested transcripts
  for (i in seq_along(transcripts)) {
    
    transcript <- transcripts[i]
    
    # get all data and lay-out features to plot current transcript
    data <- as.data.frame(cbind(t(usage[transcript,,drop=FALSE]),t(totalCount[transcript,,drop=FALSE])))
    data$group <- names(cell_to_group)[match(rownames(data),cell_to_group)]
    colnames(data) <- c("usage","totalCount","group")
    data$group <- as.factor(data$group)
    data$usage <- as.numeric(data$usage)
    data$totalCount <- as.numeric(data$totalCount)
    data$variable <- transcript # for adding Tx as a label to the facet
    padj <- format(topTable[transcript,"empirical_FDR"], digits=4)
    gene <- tx2gene[tx2gene$isoform_id == transcript,"gene_id"]
    
    # violin plot
    gg <- data %>%
      ggplot(aes(x=group,y=usage,fill=group,width=totalCount)) +
      geom_violin()  +
      geom_point(data = data, aes(x=group,y=usage, size=totalCount), position = position_jitterdodge(jitter.width = 0.7,jitter.height = 0, dodge.width = 0.9)) + scale_radius(name = "expression",range = c(0, 5)) +
      ylim(c(-0.05,1.05)) +
      ylab("proportion of usage") +
      theme_bw() +
      ggtitle(paste0(transcript, "  -  ", gene)) +
      theme(plot.title = element_text(size = 9.5,face="bold")) +
      facet_wrap(~variable, ncol=1, labeller = labeller(variable = label_facet(data$variable, padj))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8)) +
      theme(strip.text = element_text(size = 7,face="bold"))
    
    # add summarystats
    if ("model" %in% summaryStat) {
      
      model_estimates <- rowData(object)[["fitQBModels"]][transcript][[1]]@params$coefficients
      
      coefs <- do.call(rbind,coefficients)
      requested_estimates <- boot::inv.logit(coefs %*% model_estimates)
      
      data$modelled_value <- data$group
      levels(data$modelled_value) <- requested_estimates
      data$modelled_value <- as.numeric(as.character(data$modelled_value))
      
      gg <- gg +
        geom_point(data = data, aes(x=group,y=modelled_value), position = position_jitterdodge(jitter.width = 0,jitter.height = 0, dodge.width = 0.9), shape=18, size=5, colour = "gold2") 
    }
    
    if("mean"%in%summaryStat){
      gg <- gg + stat_summary(fun=mean, geom="point", position = position_dodge(width = 0.9), shape=18, size=4, colour = "cyan")
    }
    
    if("median"%in%summaryStat){
      gg <- gg + stat_summary(fun=median, geom="point", position = position_dodge(width = 0.9), shape=18, size=3, colour = "darkgreen")
    }
    
    # save current ggplot in list
    plotList[[i]] <- gg
    
  }
  return(plotList)
}

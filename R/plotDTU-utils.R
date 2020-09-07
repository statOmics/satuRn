getTotalCount <- function(countData, tx2gene){
  # get tx2gene in better format
  geneForEachTx <- tx2gene$gene_id[match(rownames(countData),tx2gene$isoform_id)]
  geneForEachTx <- as.character(geneForEachTx)
  stopifnot(class(geneForEachTx) %in% c("character", "factor"))
  stopifnot(length(geneForEachTx) == nrow(countData))

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

label_facet <- function(txID, p_adj){
  lev <- levels(as.factor(txID))
  lab <- paste0(lev, ": empFDR = ", p_adj)
  names(lab) <- lev
  return(lab)
}

visualize_DTU <- function(object, topTable, contrast, summaryStat, tx2gene, transcripts) {

  plotList <- list() # list to store ggplot objects

  countData <- assay(object)

  metaData <- as.data.frame(colData(object))
  metaData <- metaData[metaData$group %in% rownames(contrast)[as.logical(contrast)],]

  ## To calculate usages, I also need the info of transcripts that I don't even want to visualize
  genes <- unique(tx2gene[match(transcripts,tx2gene$isoform_id), "gene_id"])

  # Get all transcripts belonging  to any of the genes of interest
  transcripts_all <- c()
  for (gene in genes) {
    transcripts_all <- c(transcripts_all,as.character(tx2gene[which(tx2gene$gene_id == gene), "isoform_id"]))
  }

  # select the required cells and transcripts from countData and tx2gene
  countData <- countData[transcripts_all,metaData$sample_name]
  tx2gene <- tx2gene[tx2gene$isoform_id %in% transcripts_all,]

  # get total counts (gene-level count)
  totalCount <- getTotalCount(countData, tx2gene)

  # get usage (relative proportion of expression)
  usage <- countData/totalCount

  # subset to only the transcripts we want to plot
  usage <- usage[transcripts,,drop=FALSE]
  tx2gene <- tx2gene[tx2gene$isoform_id %in% transcripts,]
  totalCount <- totalCount[transcripts,,drop=FALSE]

  for (i in 1:length(transcripts)) {

    transcript <- transcripts[i]

    data_hlp <- cbind(metaData,t(totalCount[transcript,,drop=FALSE]))
    data_hlp <- suppressWarnings(melt(data_hlp))

    data <- cbind(metaData,t(usage[transcript,,drop=FALSE]))
    data <- suppressWarnings(melt(data))
    data$totalCount <- data_hlp$value

    padj <- format(topTable[transcript,"empirical_FDR"], digits=4)
    gene <- tx2gene[tx2gene$isoform_id == transcript,"gene_id"]

    sizeVector <- data_hlp$value
    sizeVector_scale <- rescale(sizeVector, c(0.5,5))

    # for the violin plots and general lay-out
    gg <- data %>%
      ggplot(aes(x=cluster,y=value,fill=brain_region,width=totalCount)) +
      geom_violin() +
      scale_fill_manual(values=c("royalblue4", "firebrick")) +
      ylim(c(-0.05,1.05)) +
      ylab("proportion of usage") +
      theme_bw() +
      ggtitle(gene) +
      theme(plot.title = element_text(size = 12,face="bold")) +
      facet_wrap(~variable, ncol=1, labeller = labeller(variable = label_facet(data$variable, padj))) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=6)) +
      theme(strip.text = element_text(size = 6)) +
      geom_point(data = data, aes(x=cluster,y=value), position = position_jitterdodge(jitter.width = 0.7,jitter.height = 0, dodge.width = 0.9), size=sizeVector_scale)

    if ("mean" %in% summaryStat) {
      gg  <- gg +
        stat_summary(fun=mean, geom="point", position = position_dodge(width = 0.9), shape=18, size=6, colour = "cyan")
    }

    if ("median" %in% summaryStat) {
      gg  <- gg +
        stat_summary(fun=median, geom="point", position = position_dodge(width = 0.9), shape=18, size=8, colour = "darkgreen")
    }

    if ("weightedMean" %in% summaryStat) {

      group1 <- rownames(contrast)[as.logical(contrast)][1]
      group2 <- rownames(contrast)[as.logical(contrast)][2]

      val_1 <- sum(data$value[which(metaData$group == group1)]*data$totalCount[which(metaData$group == group1)],na.rm = T)/sum(data$totalCount[which(metaData$group == group1)],na.rm = T)
      val_2 <- sum(data$value[which(metaData$group == group2)]*data$totalCount[which(metaData$group == group2)],na.rm = T)/sum(data$totalCount[which(metaData$group == group2)],na.rm = T)

      data$weighted_value <- NA
      data$weighted_value[which(metaData$group == group1)] <- val_1
      data$weighted_value[which(metaData$group == group2)] <- val_2

      gg <- gg +
        geom_point(data = data, aes(x=cluster,y=weighted_value), position = position_jitterdodge(jitter.width = 0,jitter.height = 0, dodge.width = 0.9), shape=18, size=6, colour = "purple")

    }

    if ("model" %in% summaryStat) {

      models_glm <- rowData(object)[["fitQBModels"]]

      group1 <- rownames(contrast)[as.logical(contrast)][1]
      group2 <- rownames(contrast)[as.logical(contrast)][2]

      val_1 <- boot::inv.logit(unname(models_glm[transcript][[1]]@params$coefficients[paste0("designgroup",group1)]))
      val_2 <- boot::inv.logit(unname(models_glm[transcript][[1]]@params$coefficients[paste0("designgroup",group2)]))

      data$modelled_value <- NA
      data$modelled_value[which(metaData$group == group1)] <- val_1
      data$modelled_value[which(metaData$group == group2)] <- val_2

      gg <- gg +
        geom_point(data = data, aes(x=cluster,y=modelled_value), position = position_jitterdodge(jitter.width = 0,jitter.height = 0, dodge.width = 0.9), shape=18, size=4, colour = "gold")
    }

    plotList[[i]] <- gg
  }
  return(plotList)
}

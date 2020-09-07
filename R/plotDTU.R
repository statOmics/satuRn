#' Plot function to visualize differential transcript usage (DTU)
#'
#' @description Plot function to visualize differential transcript usage (DTU)
#'
#' @param object A `SummarizedExperiment` containing the models and results
#'      of the DTU analysis as obtained by the `fitQB` and `toptable` function from
#'      this `qbDTU` package, respectively.
#'
#' @param contrast Specifies the specific contrast for which the visualization should be
#'      returned.
#'
#' @param summaryStat Which summaryStatistic for the relative usage of the transcript
#'      should be displayed. `Character` or `character vector`, must be any of following
#'      summary statistics; model (default), mean, median and/or weighted mean.
#'
#' @param transcripts A `character` or `character vector` of transcript IDs,
#'     to specify which transcripts should be visualized. Can be used together with
#'     `genes`. If not specified, `plotDTU` will check if the `genes` slot is
#'     specified.
#'
#' @param genes A single `character` or `character vector` of gene IDs,
#'     to specify the genes for which the individual transcripts should be visualized.
#'     Can be used together with `transcripts`. If not specified, `plotDTU` will check
#'     if the `transcripts` slot is specified.
#'
#' @param top.n A `numeric` value. If neither `transcripts` nor `genes` was
#'     was specified, this argument leads to the visualization of the `n` most
#'     significant DTU transcripts in the contrast. Defaults to 6 transcripts
#'
#' @examples #TODO
#'
#' @return A plot directly displayed in the current R session.
#'
#' @rdname plotDTU
#'
#' @author Jeroen Gilis
#'
#' @import ggplot2
#'
#' @export

plotDTU <- function(object, contrast, summaryStat = "model", transcripts = NULL, genes = NULL, top.n = 6){

  ## Stop if some input is not provided or not in the correct format
  stopifnot(class(object) == "SummarizedExperiment")
  stopifnot(class(contrast)[1] == "matrix")
  stopifnot(class(transcripts) %in% c("character", "NULL"))
  stopifnot(class(genes) %in% c("character", "NULL"))
  stopifnot(class(top.n) %in% c("numeric"))

  topTable <- rowData(object)[[paste0("fitQBResult_",colnames(contrast))]]
  topTable <- topTable[order(topTable$empirical_pval),]

  tx2gene <- data.frame(cbind(rowData(object)[["isoform_id"]],rowData(object)[["gene_id"]]))
  colnames(tx2gene) <- c("isoform_id","gene_id")
  tx2gene$isoform_id <- as.character(tx2gene$isoform_id)
  tx2gene$gene_id <- as.character(tx2gene$gene_id)

  # If both transcripts and genes are null
  if (is.null(transcripts) & is.null(genes)){
    transcripts <- rownames(topTable)[1:top.n]
    return(visualize_DTU(object, topTable, contrast, summaryStat, tx2gene, transcripts))
  }

  tx_tx <- tx_gene <- c()

  # If transcripts is not null
  # By stating this first, transcripts have priority over genes
  if (!is.null(transcripts)){
    absent <- NULL
    ## check if all provided transcripts are present in topTable
    absent <- transcripts[!transcripts %in% rownames(topTable)]
    if(length(absent)>0){
      warning(paste("The requested transcript", absent, "is not present in the provided topTable. "))
    }
    tx_tx <- transcripts[!transcripts %in% absent]
  }

  # If genes is not null
  if (!is.null(genes)){
    absent <- NULL
    absent <- genes[!genes %in% tx2gene$gene_id]
    if(length(absent)>0){
      warning(paste("The requested gene", absent, "is not present in the provided tx2gene dataframe. "))
    }
    genes <- genes[!genes %in% absent]
    tx_gene <- tx2gene[tx2gene$gene_id %in% genes, "isoform_id"]
  }

  transcripts <- c(tx_tx, tx_gene)

  if (length(transcripts) < 1){
    stop("None of the requested transcripts/genes could be retrieved from the provided data")
  }

  plotList <- visualize_DTU(object, topTable, contrast, summaryStat, tx2gene, transcripts)

  return(plotList)
}

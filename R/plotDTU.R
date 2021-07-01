# Computes the gene-level count in each sample or cell in the data.
# Essentially taking the sum of the counts for each transcript within a gene.
getTotalCount <- function(countData, tx2gene) {
    # get tx2gene in better format
    geneForEachTx <- tx2gene$gene_id[match(rownames(countData), 
                                            tx2gene$isoform_id)]
    geneForEachTx <- as.character(geneForEachTx)
    stopifnot(class(geneForEachTx) %in% c("character", "factor"))
    stopifnot(length(geneForEachTx) == nrow(countData))

    # adapted from DEXSeq source code
    forCycle <- split(seq_len(nrow(countData)), as.character(geneForEachTx))
    all <- lapply(forCycle, function(i) {
        sct <- countData[i, , drop = FALSE]
        rs <- t(vapply(seq_len(nrow(sct)), function(r) 
            colSums(sct[, , drop = FALSE]), numeric(ncol(countData))))
        rownames(rs) <- rownames(sct)
        rs
    })

    totalCount <- do.call(rbind, all)
    totalCount <- totalCount[rownames(countData), ]
    return(totalCount)
}

# Helper function for creating the facet header in the figure
label_facet <- function(txID, padj) {
    lev <- levels(as.factor(txID))
    lab <- paste0("satuRn empFDR = ", padj)
    names(lab) <- lev
    return(lab)
}

# Worker function of the plotDTU wrapper function
plotDTU_internal <- function(object, topTable, contrast, coefficients, groups, 
                            summaryStat, tx2gene, transcripts) {
    group <- modelled_value <- NULL
    plotList <- list()
    countData <- assay(object)

    # To calculate usages, we also need the info of transcripts that we do not 
    # want to visualize, i.e. all transcripts of all corresponding genes
    genes <- unique(tx2gene[match(transcripts, tx2gene$isoform_id), "gene_id"])
    transcripts_all <- tx2gene[tx2gene$gene_id %in% genes, "isoform_id"]

    # assign cells to a certain group (violin)
    names(groups) <- seq_along(groups)
    cell_to_group <- unlist(groups)
    names(cell_to_group) <- rep(paste0("violin", names(groups)), 
                                times = lengths(groups))

    countData <- countData[transcripts_all, unlist(groups)]
    tx2gene <- tx2gene[tx2gene$isoform_id %in% transcripts_all, ]

    # get total counts (gene-level count) and 
    # usages (relative proportion of expression)
    totalCount <- getTotalCount(countData, tx2gene)
    usage <- countData / totalCount

    # subset to only the transcripts we want to plot
    usage <- usage[transcripts, , drop = FALSE]
    totalCount <- totalCount[transcripts, , drop = FALSE]

    # plot each of the requested transcripts
    for (i in seq_along(transcripts)) {
        transcript <- transcripts[i]

        # get all data and lay-out features to plot current transcript
        data <- as.data.frame(cbind(t(usage[transcript, , drop = FALSE]), 
                                    t(totalCount[transcript, , drop = FALSE])))
        data$group <- names(cell_to_group)[match(rownames(data), cell_to_group)]
        colnames(data) <- c("usage", "totalCount", "group")
        data$group <- as.factor(data$group)
        data$usage <- as.numeric(data$usage)
        data$totalCount <- as.numeric(data$totalCount)
        data$variable <- transcript # for adding Tx as a label to the facet
        padj <- format(topTable[transcript, "empirical_FDR"], digits = 4)
        gene <- tx2gene[tx2gene$isoform_id == transcript, "gene_id"]

        # violin plot
        gg <- ggplot(data = data, 
                    aes(x = group, y = usage, fill = group, 
                        width = totalCount)) +
            geom_violin() +
            geom_point(data = data, 
                    aes(x = group, y = usage, size = totalCount), 
                    position = position_jitterdodge(jitter.width = 0.7, 
                                                    jitter.height = 0, 
                                                    dodge.width = 0.9)) +
            scale_radius(name = "expression", range = c(0, 5)) +
            ylim(c(-0.05, 1.05)) +
            ylab("Fraction of usage") +
            theme_bw() +
            ggtitle(paste0(transcript, "  -  ", gene)) +
            theme(plot.title = element_text(size = 9.5, face = "bold")) +
            facet_wrap(~variable, 
                        ncol = 1, 
                        labeller = labeller(variable = label_facet(
                            data$variable, padj))) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                                hjust = 1, size = 8)) +
            theme(strip.text = element_text(size = 7, face = "bold"))

        # add summarystats
        if ("model" %in% summaryStat) {
            model_estimates <- rowData(
                object)[["fitDTUModels"]][transcript][[1]]@params$coefficients

            coefs <- do.call(rbind, coefficients)
            requested_estimates <- boot::inv.logit(coefs %*% model_estimates)

            data$modelled_value <- data$group
            levels(data$modelled_value) <- requested_estimates
            data$modelled_value <- as.numeric(as.character(data$modelled_value))

            gg <- gg +
                geom_point(data = data, 
                        aes(x = group, y = modelled_value), 
                            position = position_jitterdodge(jitter.width = 0, 
                                                            jitter.height = 0, 
                                                            dodge.width = 0.9), 
                            shape = 18, size = 5, colour = "gold2")
        }

        if ("mean" %in% summaryStat) {
            gg <- gg + stat_summary(fun = mean, geom = "point", 
                                    position = position_dodge(width = 0.9), 
                                    shape = 18, size = 4, colour = "cyan")
        }

        if ("median" %in% summaryStat) {
            gg <- gg + stat_summary(fun = median, geom = "point", 
                                    position = position_dodge(width = 0.9), 
                                    shape = 18, size = 3, colour = "darkgreen")
        }

        # save current ggplot in list
        plotList[[i]] <- gg
    }
    return(plotList)
}

#' Plot function to visualize differential transcript usage (DTU)
#'
#' @description Plot function to visualize differential transcript usage (DTU)
#'
#' @param object A `SummarizedExperiment` containing the models and results
#'     of the DTU analysis as obtained by the `fitDTU` and `testDTU` function 
#'     from this `satuRn` package, respectively.
#'
#' @param contrast Specifies the specific contrast for which the visualization 
#'     should be returned. This should be one of the column names of 
#'     the contrast matrix that was provided to the `testDTU` function.
#'
#' @param groups A `list` containing two character vectors.
#'     Each character vector contains the names (sample names or cell names) 
#'     of the respective groups in the target contrast.
#'
#' @param coefficients A `list` containing two numeric vectors. 
#'     Each numeric vector specifies the model coefficient of the corresponding 
#'     groups in the selected contrast.
#'
#' @param summaryStat Which summary statistic for the relative usage of 
#'     the transcript should be displayed. `Character` or `character vector`, 
#'     must be any of following summary statistics; model (default), 
#'     mean or median.
#'
#' @param transcripts A `character` or `character vector` of transcript IDs,
#'     to specify which transcripts should be visualized. Can be used together 
#'     with `genes`. If not specified, `plotDTU` will check if the `genes` slot 
#'     is specified.
#'
#' @param genes A single `character` or `character vector` of gene IDs,
#'     to specify the genes for which the individual transcripts should be 
#'     visualized. Can be used together with `transcripts`. If not specified, 
#'     `plotDTU` will check if the `transcripts` slot is specified.
#'
#' @param top.n A `numeric` value. If neither `transcripts` nor `genes` was
#'     was specified, this argument leads to the visualization of the `n` most
#'     significant DTU transcripts in the contrast. Defaults to 6 transcripts.
#'
#' @examples
#' data(sumExp_example, package = "satuRn")
#' library(SummarizedExperiment)
#' sumExp <- fitDTU(
#'     object = sumExp_example,
#'     formula = ~ 0 + group,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     verbose = TRUE
#' )
#' group <- as.factor(colData(sumExp)$group)
#' design <- model.matrix(~ 0 + group)
#' colnames(design) <- levels(group)
#' L <- matrix(0, ncol = 2, nrow = ncol(design))
#' rownames(L) <- colnames(design)
#' colnames(L) <- c("Contrast1", "Contrast2")
#' L[c("VISp.L5_IT_VISp_Hsd11b1_Endou", "ALM.L5_IT_ALM_Tnc"), 1] <- c(1, -1)
#' L[c("VISp.L5_IT_VISp_Hsd11b1_Endou", 
#'     "ALM.L5_IT_ALM_Tmem163_Dmrtb1"), 2] <- c(1, -1)
#'
#' sumExp <- testDTU(object = sumExp, 
#'     contrasts = L, 
#'     diagplot1 = FALSE,
#'     diagplot2 = FALSE,
#'     sort = FALSE)
#'
#' group1 <- rownames(colData(sumExp))[colData(sumExp)$group == 
#'                                     "VISp.L5_IT_VISp_Hsd11b1_Endou"]
#' group2 <- rownames(colData(sumExp))[colData(sumExp)$group == 
#'                                     "ALM.L5_IT_ALM_Tnc"]
#'
#' plots <- plotDTU(
#'     object = sumExp,
#'     contrast = "Contrast1",
#'     groups = list(group1, group2),
#'     coefficients = list(c(0, 0, 1), c(0, 1, 0)),
#'     summaryStat = "model",
#'     transcripts = c("ENSMUST00000165123", 
#'                     "ENSMUST00000165721", 
#'                     "ENSMUST00000005067"),
#'     genes = NULL,
#'     top.n = 6
#' )
#' @return A ggplot object that can be directly displayed in the current R 
#' session or stored in a list.
#'
#' @rdname plotDTU
#'
#' @author Jeroen Gilis
#'
#' @import ggplot2
#' @importFrom boot inv.logit
#' @importFrom Matrix Matrix
#' @importFrom stats median
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#'
#' @export

# Wrapper function, sanity checks and 
# getting all transcripts + the requested contrast in place
plotDTU <- function(object, contrast, groups, coefficients, 
                    summaryStat = "model", transcripts = NULL, genes = NULL, 
                    top.n = 6) {

    ## Stop if some input is not provided or not in the correct format
    stopifnot(is(object, "SummarizedExperiment"))
    stopifnot(class(transcripts) %in% c("character", "NULL"))
    stopifnot(class(genes) %in% c("character", "NULL"))
    stopifnot(class(top.n) %in% c("numeric"))
    if (is.null(rowData(object)[["fitDTUModels"]])) {
        stop("fitDTUModels is empty. Did you run fitDTU first?")
    }
    if (!any(grepl("fitDTUResult", names(rowData(object))))) {
        stop("fitDTUResult is empty. Did you run testDTU first?")
    }
    
    # select the requested contrast
    topTable <- rowData(object)[[paste0("fitDTUResult_", contrast)]] 
    topTable <- topTable[order(topTable$empirical_pval), ]

    tx2gene <- data.frame(cbind(rowData(object)[["isoform_id"]], 
                                rowData(object)[["gene_id"]]))
    colnames(tx2gene) <- c("isoform_id", "gene_id")
    tx2gene$isoform_id <- as.character(tx2gene$isoform_id)
    tx2gene$gene_id <- as.character(tx2gene$gene_id)

    # If both transcripts and genes are null, 
    # plot the top n DTU transcripts (default = 6)
    if (is.null(transcripts) & is.null(genes)) {
        transcripts <- rownames(topTable)[seq_len(top.n)]
        return(plotDTU_internal(object, topTable, contrast, coefficients, 
                                groups, summaryStat, tx2gene, transcripts))
    }

    tx_tx <- tx_gene <- c()

    # If transcripts is not null
    if (!is.null(transcripts)) {
        absent <- NULL
        ## check if all provided transcripts are present in topTable
        absent <- transcripts[!transcripts %in% rownames(topTable)]
        if (length(absent) > 0) {
            warning("The requested transcript", absent, 
                    "is not present in the provided topTable. ")
        }
        tx_tx <- transcripts[!transcripts %in% absent]
    }

    # If genes is not null
    if (!is.null(genes)) {
        absent <- NULL
        absent <- genes[!genes %in% tx2gene$gene_id]
        if (length(absent) > 0) {
            warning("The requested gene", absent, 
                    "is not present in the provided tx2gene dataframe. ")
        }
        genes <- genes[!genes %in% absent]
        tx_gene <- tx2gene[tx2gene$gene_id %in% genes, "isoform_id"]
    }

    transcripts <- c(tx_tx, tx_gene)

    if (length(transcripts) < 1) {
        stop("None of the requested transcripts/genes could be retrieved 
            from the provided data")
    }

    # got to the internal visualization function
    plotList <- plotDTU_internal(object, topTable, contrast, coefficients, 
                                groups, summaryStat, tx2gene, transcripts)
    return(plotList)
}

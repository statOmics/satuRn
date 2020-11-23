# Compute the dispersion parameter of the quasi-binomial GLM
.calcDispersion <- function(model, type) {
    model$dispersion <- NA
    if (type != "fitError") {
        df.r <- model$df.residual
        if (df.r > 0) {
            model$dispersion <- sum((model$weights * model$residuals^2)[model$weights > 0]) / df.r
        }
    }
    return(model)
}

# Compute the unscaled variance-covariance matrix of the quasi-binomial GLM
.calcVcovUnscaled <- function(model, type) {
    if (!type == "fitError") {
        p1 <- 1L:model$rank
        p <- length(model$coef)
        out <- matrix(NA, p, p)
        out[p1, p1] <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
        colnames(out) <- rownames(out) <- names(model$coefficients)
        model$vcovUnsc <- out
    } else {
        model$vcovUnsc <- NA
    }
    return(model)
}

# Computes the "other count" for all in each sample or cell in the data.
# Essentially taking the sum of the counts for each of the transcripts within a gene except the target transcript.
# Can be thought as the difference between the gene-level count and current transcript-level count.
.getOtherCount <- function(countData, tx2gene) {
    # get tx2gene in better format
    geneForEachTx <- tx2gene$gene_id[match(rownames(countData), tx2gene$isoform_id)]
    geneForEachTx <- as.character(geneForEachTx)
    stopifnot(class(geneForEachTx) %in% c("character", "factor"))
    stopifnot(length(geneForEachTx) == nrow(countData))
    
    forCycle <- split(1:nrow(countData), as.character(geneForEachTx))
    all <- lapply(forCycle, function(i) {
        sct <- countData[i, , drop = FALSE]
        rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[-r, , drop = FALSE])))
        rownames(rs) <- rownames(sct)
        rs
    })
    
    otherCount <- do.call(rbind, all)
    otherCount <- otherCount[rownames(countData), ]
    return(otherCount)
}

# Worker function that fits quasi-binomial models, wrapped inside the fitDTU function
.fitDTU_internal <- function(countData, tx2gene, design, parallel, BPPARAM, verbose) {
    if (parallel) {
        BiocParallel::register(BPPARAM)
        if (verbose) {
            # update progress bar 40 times
            BPPARAM$tasks <- as.integer(40)
            # show progress bar
            BPPARAM$progressbar <- TRUE
        }
    }
    
    stopifnot(class(countData)[1] %in% c("matrix", "data.frame"))
    countData <- as.matrix(countData)
    
    # Get the "other" counts, i.e. the counts for all other transcripts belonging to
    # the same gene as the current transcript
    otherCount <- satuRn:::.getOtherCount(countData, tx2gene)
    stopifnot(all(rownames(countData) %in% rownames(otherCount)))
    
    # The actual fit function
    fitQuasiLogistic <- function(countData, otherCount, design) {
        countsAll <- cbind(countData, otherCount)
        drop <- rowSums(countsAll) == 0 ## gene count is zero in this cell
        countsAll <- countsAll + 1
        countsAll[drop, ] <- NA
        
        model <- try(glm(countsAll ~ -1 + design, family = "quasibinomial"))
        
        if (class(model)[1] == "try-error") {
            model <- list()
            type <- "fitError"
            class(model) <- "list"
        } else {
            type <- "glm"
            class(model) <- "list"
        }
        
        model <- satuRn:::.calcDispersion(model, type) ## calculate disp slot
        model <- satuRn:::.calcVcovUnscaled(model, type) ## calculate vcov slot
        
        model <- model[c("coefficients", "df.residual", "dispersion", "vcovUnsc")]
        
        
        .out <- .StatModel(
            type = type,
            params = model,
            varPosterior = as.numeric(NA),
            dfPosterior = as.numeric(NA)
        )
        return(.out)
    }
    
    # Fit the models
    if (parallel) {
        models <- BiocParallel::bplapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData = countData[i, ], otherCount = otherCount[i, ], design = design), BPPARAM = BPPARAM)
    } else {
        if (verbose) {
            models <- pbapply::pblapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData = countData[i, ], otherCount = otherCount[i, ], design = design))
        } else {
            models <- lapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData = countData[i, ], otherCount = otherCount[i, ], design = design))
        }
    }
    
    # retain transcript names
    names(models) <- rownames(countData)
    
    # Squeeze a set of sample variances together by computing empirical Bayes posterior means
    hlp <- limma::squeezeVar(
        var = unlist(sapply(models, satuRn:::getDispersion)),
        df = unlist(sapply(models, satuRn:::getDF)),
        robust = FALSE
    )
    
    # put variance and degrees of freedom in appropriate slots
    for (i in 1:length(models)) {
        mydf <- hlp$df.prior + getDF(models[[i]])
        models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
        models[[i]]@dfPosterior <- as.numeric(mydf)
    }
    
    # return object of class StatModel
    return(models)
}

#' Function to fit quasi-binomial models
#'
#' @description Parameter estimation of quasi-binomial models.
#'
#' @param object A `SummarizedExperiment` instance
#'
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting.
#'
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{bpparam} in \code{BiocParallel} package for details.
#'
#' @param verbose Logical, should progress be printed?
#'
#' @examples
#' data(sumExp_vignette, package = "satuRn")
#' data(Tasic_metadata_vignette, package = "satuRn")
#' Tasic_metadata_vignette$group <- paste(Tasic_metadata_vignette$brain_region, Tasic_metadata_vignette$cluster, sep = ".")
#' sumExp <- fitDTU(
#'    object = sumExp_vignette,
#'    parallel = TRUE,
#'    BPPARAM = BiocParallel::bpparam(),
#'    verbose = TRUE)
#' 
#' @return An updated `SummarizedExperiment` instance. The instance now includes
#' a new list of models ("fitDTUModels") in its rowData slot, which can be accessed
#' by rowData(object)[["fitDTUModels"]].
#'
#' @rdname fitDTU
#'
#' @author Jeroen Gilis
#'
#' @import SummarizedExperiment
#' @importFrom limma squeezeVar
#' @importFrom BiocParallel register
#' @importFrom BiocParallel bplapply
#' @importFrom pbapply pblapply
#'
#' @export

# Wrapper function
setMethod(
    f = "fitDTU",
    signature = "SummarizedExperiment",
    function(object,
    parallel = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE) {
        if (ncol(colData(object)) == 0) stop("error: colData is empty")
        design <- model.matrix(object@metadata$formula, colData(object))

        rowData(object)[["fitDTUModels"]] <- satuRn:::.fitDTU_internal(
            countData = assay(object),
            tx2gene = rowData(object)[, c("isoform_id", "gene_id")],
            design = design,
            parallel = parallel,
            BPPARAM = BPPARAM,
            verbose = verbose
        )

        return(object)
    }
)

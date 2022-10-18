# Fit a quasibinomial regression model, given counts for all features within a
# gene and a design matrix
.fitQuasiLogistic <- function(mat, design) {
  
  if(nrow(mat) == 1){
    .out <- .StatModel(
      type = "lonelyTranscript",
      params = list(coefficients = NA, 
                    df.residual = NA,
                    dispersion = NA,
                    vcovUnsc = NA),
      varPosterior = as.numeric(NA),
      dfPosterior = as.numeric(NA)
    )
    return(.out)
  }
  
  mat <- as.matrix(mat)
  totalCount <- matrix(data = rep(colSums(mat),nrow(mat)),
                       nrow = nrow(mat),
                       byrow = TRUE)
  otherCount <- totalCount-mat
  
  models_gene <- lapply(seq_len(nrow(mat)), function(i){
    countsAll <- cbind(mat[i,], otherCount[i,])
    drop <- rowSums(countsAll) == 0 ## gene count is zero in this cell
    countsAll <- countsAll + 1
    countsAll[drop, ] <- NA
    
    model <- try(glm(countsAll ~ -1 + design, family = "quasibinomial"))
    
    # if the quasibinomial model could not be estimated, return empty model
    if (class(model)[1] == "try-error") {
      .out <- .StatModel(
        type = "fitError",
        params = list(coefficients = NA, 
                      df.residual = NA,
                      dispersion = NA,
                      vcovUnsc = NA),
        varPosterior = as.numeric(NA),
        dfPosterior = as.numeric(NA)
      )
      return(.out)
    } else {
      type <- "glm"
      class(model) <- "list"
      model <- .calcDispersion(model, type) ## calculate disp slot
      model <- .calcVcovUnscaled(model, type) ## calculate vcov slot
      
      model <- model[c("coefficients", 
                       "df.residual", 
                       "dispersion", 
                       "vcovUnsc")]
      
      .out <- .StatModel(
        type = type,
        params = model,
        varPosterior = as.numeric(NA),
        dfPosterior = as.numeric(NA)
      )
      return(.out)
    }
  })
}

# Compute the dispersion parameter of the quasi-binomial GLM
.calcDispersion <- function(model, type) {
    model$dispersion <- NA
    df.r <- model$df.residual
    if (df.r > 0) {
        model$dispersion <- sum(
            (model$weights * model$residuals^2)[model$weights > 0]) / df.r
    }
    return(model)
}

# Compute the unscaled variance-covariance matrix of the quasi-binomial GLM
.calcVcovUnscaled <- function(model, type) {
    p1 <- 1L:model$rank
    p <- length(model$coef)
    out <- matrix(NA, p, p)
    out[!is.na(model$coefficients), !is.na(model$coefficients)] <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    colnames(out) <- rownames(out) <- names(model$coefficients)
    model$vcovUnsc <- out
    return(model)
}

# Worker function that fits quasi-binomial models, 
# wrapped inside the fitDTU function
.fitDTU_internal <- function(countData, 
                             tx2gene, 
                             design, 
                             parallel, 
                             BPPARAM, 
                             verbose) {
  if (parallel) {
    BiocParallel::register(BPPARAM)
    if (verbose) {
      # update progress bar 40 times
      BPPARAM$tasks <- as.integer(40)
      # show progress bar
      BPPARAM$progressbar <- TRUE
    }
  }
  
  stopifnot(class(countData)[1] %in% c("matrix", "data.frame", 
                                       "dgCMatrix", "DelayedMatrix"))
  
  # Check if lonely transcripts (only transcript in its corresponding gene)
  # have already been removed, if not throw message
  lonely_tx <- !tx2gene$gene_id %in% tx2gene$gene_id[duplicated(tx2gene$gene_id)]
  if(any(lonely_tx)){
    message("Genes with only one type of transcript/isoform detected. 
            Results of such transcripts will be set to NA and flagged as
            lonelyTranscript.")
  }
  
  geneForEachTx <- tx2gene$gene_id[match(rownames(countData), 
                                         tx2gene$isoform_id)]
  geneForEachTx <- as.character(geneForEachTx)
  stopifnot(length(geneForEachTx) == nrow(countData))
  
  # split (sparse) matrix in a per-gene list of (sparse) matrices
  colnames(countData) <- NULL # to avoid having colnames in each sub-matrix
  matList <- split.data.frame(countData, geneForEachTx)
  
  # Fit the models
  if (parallel) {
    models <- BiocParallel::bplapply(matList, function(mat){
      models_gene <- .fitQuasiLogistic(mat = mat, design = design)
      models_gene <- unlist(models_gene)
      names(models_gene) <- rownames(mat)
      return(models_gene)
    }, BPPARAM = BPPARAM)
  } else if(verbose) {
    models <- pbapply::pblapply(matList, function(mat){
      models_gene <- .fitQuasiLogistic(mat = mat, design = design)
      models_gene <- unlist(models_gene)
      names(models_gene) <- rownames(mat)
      return(models_gene)
    })
  } else {
    models <- lapply(matList, function(mat){
      models_gene <- .fitQuasiLogistic(mat = mat, design = design)
      names(models_gene) <- rownames(mat)
      return(models_gene)
    })
  }
  
  names(models) <- NULL
  models <- unlist(models)
  models <- models[match(rownames(countData), names(models))]
  
  # Squeeze a set of sample variances together 
  # by computing empirical Bayes posterior means
  hlp <- limma::squeezeVar(
    var = unlist(vapply(models, getDispersion, numeric(1))),
    df = unlist(vapply(models, getDF, numeric(1))),
    robust = FALSE
  )
  
  # put variance and degrees of freedom in appropriate slots
  for (i in seq_along(models)) {
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
#' @param object A `SummarizedExperiment` instance generated with the
#' SummarizedExperiment function of the SummarizedExperiment package.
#' In the assay slot, provide the transcript-level expression counts as an
#' ordinary `matrix`, `DataFrame`, a `sparseMatrix` or a `DelayedMatrix`.
#' The `rowData` slot must be a `DataFrame` object describing the rows,
#' which must contain a column `isoform_id` with the row names of
#' the expression matrix and a column `gene_id` with the corresponding gene
#' identifiers of each transcript. `colData` is a `DataFrame` describing
#' the samples or cells in the experiment. Finally, specify the
#' experimental design as a formula in the metadata slot. This formula must
#' be based on the colData. See the documentation examples and the vignette
#' for more details.
#'
#' @param formula Model formula. The model is built based on the
#' covariates in the data object.
#'
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting procedure.
#'
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#' back-end to be used for computations. See
#' \code{bpparam} in \code{BiocParallel} package for details.
#'
#' @param verbose Logical, should progress be printed?
#'
#' @examples
#' data(sumExp_example, package = "satuRn")
#' sumExp <- fitDTU(
#'     object = sumExp_example,
#'     formula = ~ 0 + group,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     verbose = TRUE
#' )
#' @return An updated `SummarizedExperiment` instance. The instance now includes
#'     a new list of models ("fitDTUModels") in its rowData slot,
#'     which can be accessed by rowData(object)[["fitDTUModels"]].
#'
#' @rdname fitDTU
#'
#' @author Jeroen Gilis
#'
#' @import SummarizedExperiment
#' @importFrom limma squeezeVar
#' @importFrom stats model.matrix
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
    formula,
    parallel = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE) {
        if (ncol(SummarizedExperiment::colData(object)) == 0) 
            stop("colData is empty")

        design <- model.matrix(formula, SummarizedExperiment::colData(object))

        if (!"gene_id" %in% colnames(rowData(object)) | 
            !"isoform_id" %in% colnames(rowData(object))) {
            stop("rowData does not contain columns gene_id and isoform_id")
        }
        if (!all(rownames(object) == rowData(object)[, "isoform_id"])) {
            stop("not all row names of the expression matrix match 
                the isoform_id column of the object's rowData")
        }

        rowData(object)[["fitDTUModels"]] <- .fitDTU_internal(
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

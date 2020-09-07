#' Function to fit quasi-binomial models
#'
#' @description Parameter estimation of quasi-binomial models.
#'
#' @param object `SummarizedExperiment` instance
#'
#' @param modelColumnName `character` to indicate the variable name that is used
#'        to store the fitQB models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the Features instance. Default is "fitQBModels".
#'
#' @examples #TODO
#'
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname fitQB
#'
#' @author Jeroen Gilis
#'
#' @import SummarizedExperiment
#' @importFrom speedglm speedglm
#' @importFrom limma squeezeVar
#' @importFrom BiocParallel register
#' @importFrom BiocParallel bplapply
#' @importFrom pbapply pblapply
#'
#' @export

setMethod("fitQB","SummarizedExperiment",
      function(object,
              speed = FALSE,
              parallel = FALSE,
              BPPARAM = BiocParallel::bpparam(),
              verbose = TRUE){
          if (ncol(colData(object))==0) stop("error: colData is empty")
          design <- model.matrix(object@metadata$formula, colData(object))

          rowData(object)[["fitQBModels"]] <- fitQB_internal(countData = assay(object),
                                                      tx2gene = rowData(object)[,c("isoform_id", "gene_id")],
                                                      design = design,
                                                      speed = speed,
                                                      parallel = parallel,
                                                      BPPARAM = BPPARAM,
                                                      verbose = verbose)

          return(object)
})











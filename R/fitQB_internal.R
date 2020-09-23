#' Function to fit quasi-binomial models
#'
#' @description Parameter estimation of quasi-binomial models.
#'
#' @param countData A `DataFrame` or `matrix` of transcript expression values.
#'     The transcripts are along the rows and samples/cells along the columns.
#'
#' @param tx2gene A `DataFrame` or `matrix` that links each transcript in
#'     the countData object to its corresponding gene ID.
#'
#' @param design A `DataFrame` with information on the design. It has
#'     the same number of rows as the number of columns (samples/cells) of
#'     `countData`.
#'
#' @param speed `boolean(1)` Logical, defaults to FALSE. Set to TRUE if
#'     the quasi-binomial models should be fitted with the \code{speedglm} package
#'     rather than the classical glm procedure, which depending on the dataset might
#'     reduce computational time and memory usage. See the \code{speedglm} package
#'     for details.
#'
#' @param parallel `boolean(1)` Logical, defaults to FALSE. Set to TRUE if to
#'     parallellize the fitting procedure.
#'
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'     back-end to be used for computations. See
#'     \code{bpparam} in \code{BiocParallel} package for details.
#'
#' @param verbose `boolean(1)` Logical, defaults to TRUE. Used to indicate if a
#'      progress bar should be printed during the fitting procedure.
#'
#' @examples #TODO
#'
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname fitQB_internal
#'
#' @author Jeroen Gilis
#'
#' @importFrom speedglm speedglm
#' @importFrom limma squeezeVar
#' @importFrom BiocParallel register
#' @importFrom BiocParallel bplapply
#' @importFrom pbapply pblapply
#'
#' @export

fitQB_internal <- function(countData, tx2gene, design, speed, parallel, BPPARAM, verbose){

  if (parallel) {
    BiocParallel::register(BPPARAM)
    if (verbose) {
      # update progress bar 40 times
      BPPARAM$tasks = as.integer(40)
      # show progress bar
      BPPARAM$progressbar = TRUE
    }
  }

  stopifnot(class(countData)[1] %in% c("matrix", "data.frame"))
  countData <- as.matrix(countData)

  # Get the "other" counts, i.e. the counts for all other transcripts belonging to
  # the same gene as the current transcript
  otherCount <- getOtherCount(countData, tx2gene)
  stopifnot(all(rownames(countData) %in% rownames(otherCount)))

  # The actual fit function
  fitQuasiLogistic <- function(countData, otherCount, design, speed = speed){

    countsAll <- cbind(countData, otherCount)
    drop <- rowSums(countsAll) == 0 ## gene count is zero in this cell
    countsAll <- countsAll + 1
    countsAll[drop,] <- NA

    if (speed == FALSE) {
      model <- try(glm(countsAll ~ -1+design, family="quasibinomial"))
    } else {
      model <- try(speedglm(countsAll ~ -1+design, family=quasibinomial()))
    }

    if (class(model)[1]=="try-error"){
      model <- list()
      type <- "fitError"
      class(model) <- "list"
    } else {
      type <- "glm"
      class(model) <- "list"
    }

    if (speed == FALSE){
      model <- calcDispersion(model,type) ## calculate disp slot
      model <- calcVcovUnscaled(model,type) ## calculate vcov slot

      model <- model[c("coefficients","df.residual", "dispersion", "vcovUnsc")]
    } else { # if speed == TRUE
      if (type != "fitError") {
        model$vcovUnsc <- solve(model$XTX)
      } else {
        model$vcovUnsc <- NA
      }

      model <- model[c("coefficients","df", "dispersion", "vcovUnsc")]
      names(model) <- c("coefficients","df.residual", "dispersion", "vcovUnsc")
    }

    .out  <- .StatModel(type = type,
                        params = model,
                        varPosterior = as.numeric(NA),
                        dfPosterior = as.numeric(NA))
    return(.out)
  }

  # Fit the models
  if (parallel) {
    models <- BiocParallel::bplapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData=countData[i,], otherCount=otherCount[i,], design = design, speed = speed), BPPARAM = BPPARAM
    )
  } else {
    if (verbose) {
      models <- pbapply::pblapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData=countData[i,], otherCount=otherCount[i,], design = design, speed = speed)
      )
    } else {
      models <- lapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData=countData[i,], otherCount=otherCount[i,], design = design, speed = speed)
      )
    }
  }

  # retain transcript names
  names(models) <- rownames(countData)

  # Squeeze a set of sample variances together by computing empirical Bayes posterior means
  hlp <- limma::squeezeVar(var = sapply(models, getDispersion),
                           df = sapply(models, getDF),
                           robust = FALSE)

  # put variance and degrees of freedom in appropriate slots
  for (i in 1:length(models)) {
    mydf <- hlp$df.prior + getDF(models[[i]])
    models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
    models[[i]]@dfPosterior <- as.numeric(mydf)
  }

  # return object of class StatModel
  return(models)
}

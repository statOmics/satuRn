#' Accessor functions for StatModel class
#'
#' @description Accessor functions for StatModel class
#'              \describe{
#'              \item{getModel(object)}{to get model}
#'              \item{getDF(object)}{to get the residual degrees of freedom of the model}
#'              \item{getDfPosterior(object)}{to get the degrees of freedom of
#'                the empirical Bayes variance estimator}
#'              \item{getDispersion(object)}{to get the dispersion estimate of the model}
#'              \item{getCoef(object)}{to get the parameter estimates of the mean model}
#'              }
#'
#' @rdname statModelAccessors
#' @aliases statModelAccessors getModel getDF getDfPosterior getDispersion getCoef
#'
#' @param object `StatModel` object

setMethod("getModel",
    signature = "StatModel",
    definition = function(object) object@params
)

setMethod("getDF",
    signature = "StatModel",
    definition = function(object) object@params$df.residual
)

setMethod("getDfPosterior",
    signature = "StatModel",
    definition = function(object) object@dfPosterior
)

setMethod("getDispersion",
    signature = "StatModel",
    definition = function(object) object@params$dispersion
)

setMethod("getCoef",
    signature = "StatModel",
    definition = function(object) object@params$coefficients
)

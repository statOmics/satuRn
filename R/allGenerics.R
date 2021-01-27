#' @export
setGeneric("getModel", function(object) standardGeneric("getModel"))

#' @export
setGeneric("getDF", function(object) standardGeneric("getDF"))

#' @export
setGeneric("getDfPosterior", function(object) standardGeneric("getDfPosterior"))

#' @export
setGeneric("getDispersion", function(object) standardGeneric("getDispersion"))

#' @export
setGeneric("getCoef", function(object) standardGeneric("getCoef"))

#' @export
#' @name fitDTU
#' @title fitDTU
#' @param ... parameters including:
setGeneric("fitDTU", function(object, ...) standardGeneric("fitDTU"))

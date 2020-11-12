#' classes for satuRn
#'
#' @title The StatModel class for satuRn
#' @slot type The type of the used model
#' @slot params A list containing information of the used model
#' @slot varPosterior A vector of posterior variance
#' @slot dfPosterior A vector of posterior degrees of freedom
#' @rdname satuRn-framework
#' @author Jeroen Gilis
#' @export
.StatModel <- setClass("StatModel",
                       slots = c(type = "character",
                                 params = "list",
                                 varPosterior = "numeric",
                                 dfPosterior = "numeric"))

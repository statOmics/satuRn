##' @export
StatModel <- function(type="fitError",params=list(),varPosterior=is.numeric(NA),dfPosterior=is.numeric(NA)) {
  out<-new("StatModel")
  out@type=type
  out@params=params
  out@varPosterior=varPosterior
  out@dfPosterior=dfPosterior
  return(out)
}

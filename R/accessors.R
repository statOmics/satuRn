setMethod("getModel",signature="StatModel",definition=function(object)
  object@params)

setMethod("getDF",signature="StatModel",definition=function(object)
  object@params$df.residual)

setMethod("getDfPosterior",signature="StatModel",definition=function(object)
  object@dfPosterior)

setMethod("getDispersion",signature="StatModel",definition=function(object)
  object@params$dispersion)

setMethod("getCoef",signature="StatModel",definition=function(object)
  object@params$coefficients)

getEstimates <- function(object,contrast){
  coef <- getCoef(object)
  if (is.null(coef)){
    coef <- rep(NA, times=nrow(contrast))
  }
  return(contrast%*%coef)
}

varContrast <- function(object,contrast){
  if (object@type!="fitError"){
    if (nrow(object@params$vcovUnsc) == length(contrast)){
      vcovTmp <- object@params$vcovUnsc*object@varPosterior
      return(diag(t(contrast)%*%vcovTmp%*%contrast))
    }
  }
  return(NA)
}

## code based on source code of locFDR package https://CRAN.R-project.org/package=locfdr
p.adjust_empirical_hlp <- function(zz){

  N <- length(zz)
  b <- 4.3 * exp(-0.26*log(N,10))
  med <- median(zz)
  sc <- diff(quantile(zz)[c(2,4)])/(2*qnorm(.75))
  mlests <- locfdr:::locmle(zz, xlim=c(med, b*sc))

  nulltype = 1
  lo <- min(zz)
  up <- max(zz)
  bre = 120
  breaks <- seq(lo, up, length = bre)
  zzz <- pmax(pmin(zz, up), lo)
  zh <- hist(zzz, breaks = breaks, plot = F)
  x <- (breaks[-1] + breaks[ - length(breaks)])/2
  sw <- 0
  X <- cbind(1, poly(x, df = 7))
  zh <- hist(zzz, breaks = breaks, plot = F)
  y <- zh$counts
  f <- glm(y ~ poly(x, df = 7), poisson)$fit

  Cov.in = list(x=x, X=X, f=f, sw=sw)
  ml.out = locfdr:::locmle(zz, xlim = c(mlests[1], b * mlests[2]),
                           d=mlests[1], s=mlests[2], Cov.in=Cov.in)
  mlests = ml.out$mle

  return(mlests)

}

## code based on source code of locFDR package https://CRAN.R-project.org/package=locfdr
p.adjust_empirical <- function(pvalues,tvalues,plot=FALSE){

  zvalues <- qnorm(pvalues/2)*sign(tvalues)

  zvalues_mid <- zvalues[abs(zvalues) < 10]
  zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)] ## to avoid locFDR (numerical) fit errors if abs(z) is extremely large

  mlests <- p.adjust_empirical_hlp(zvalues_mid)

  zval_empirical <- (zvalues-mlests[1])/mlests[2]
  pval_empirical <- 2*pnorm(-abs(zval_empirical), mean=0, sd=1)

  if (plot){
    zval_empirical <- zval_empirical[!is.na(zval_empirical)]
    lo <- min(zval_empirical)
    up <- max(zval_empirical)

    lo <- min(lo, -1*up) # to center the figure (is this a good idea?)
    up <- max(up, -1*lo)

    bre = 120
    breaks <- seq(lo, up, length = bre)
    zzz <- pmax(pmin(zval_empirical, up), lo)
    zh <- hist(zzz, breaks = breaks, plot = F)
    yall <- zh$counts
    K <- length(yall)

    hist(zzz, breaks = breaks, xlab = "z-scores", main = "Empirical distribution of z-scores", freq = F)
    xfit <- seq(min(zzz), max(zzz), length = 4000)
    yfit <- dnorm(xfit/mlests[3], mean = 0, sd = 1)
    lines(xfit, yfit, col = "darkgreen", lwd = 2)
  }

  FDR <- p.adjust(pval_empirical, method = "BH")
  newList <- list("pval" = pval_empirical, "FDR" = FDR)
  return(newList)
}

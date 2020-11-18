# Compute usage estimates (log-odds) in a contrast of interest
getEstimates <- function(object, contrast) {
    coef <- getCoef(object)
    if (is.null(coef)) {
        coef <- rep(NA, times = length(contrast))
    }
    return(contrast %*% coef)
}

# Compute the variance oon the usage estimates in a contrast of interest
varContrast <- function(object, contrast) {
    if (object@type != "fitError") {
        if (nrow(object@params$vcovUnsc) == length(contrast)) {
            vcovTmp <- object@params$vcovUnsc * object@varPosterior
            return(diag(t(contrast) %*% vcovTmp %*% contrast))
        }
    }
    return(NA)
}

# Compute p-values under an empirical null distribution
# The code is based on the source code of locFDR package https://CRAN.R-project.org/package=locfdr
p.adjust_empirical <- function(pvalues, tvalues, plot = FALSE) {
    zvalues <- qnorm(pvalues / 2) * sign(tvalues)

    zvalues_mid <- zvalues[abs(zvalues) < 10]
    zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)] # to avoid locFDR (numerical) fit errors if abs(z) is extremely large

    #### start code from locFDR package
    N <- length(zvalues_mid)
    b <- 4.3 * exp(-0.26 * log(N, 10))
    med <- median(zvalues_mid)
    sc <- diff(quantile(zvalues_mid)[c(2, 4)]) / (2 * qnorm(.75))
    mlests <- locfdr:::locmle(zvalues_mid, xlim = c(med, b * sc)) # initial MLEs of empirical null parameters

    lo <- min(zvalues_mid)
    up <- max(zvalues_mid)
    bre <- 120
    breaks <- seq(lo, up, length = bre)
    zzz <- pmax(pmin(zvalues_mid, up), lo)
    zh <- hist(zzz, breaks = breaks, plot = F)
    x <- (breaks[-1] + breaks[-length(breaks)]) / 2
    sw <- 0
    X <- cbind(1, poly(x, df = 7))
    zh <- hist(zzz, breaks = breaks, plot = F)
    y <- zh$counts
    f <- glm(y ~ poly(x, df = 7), poisson)$fit

    Cov.in <- list(x = x, X = X, f = f, sw = sw)
    ml.out <- locfdr:::locmle(zvalues_mid,
        xlim = c(mlests[1], b * mlests[2]),
        d = mlests[1], s = mlests[2], Cov.in = Cov.in
    ) # updated MLEs of empirical null parameters
    mlests <- ml.out$mle
    #### end code from locFDR package

    zval_empirical <- (zvalues - mlests[1]) / mlests[2]
    pval_empirical <- 2 * pnorm(-abs(zval_empirical), mean = 0, sd = 1) # p-values under empirical null

    if (plot) {
        zval_empirical <- zval_empirical[!is.na(zval_empirical)]
        lo <- min(zval_empirical)
        up <- max(zval_empirical)

        lo <- min(lo, -1 * up) # to center the figure (is this a good idea?)
        up <- max(up, -1 * lo)

        bre <- 120
        breaks <- seq(lo, up, length = bre)
        zzz <- pmax(pmin(zval_empirical, up), lo)
        zh <- hist(zzz, breaks = breaks, plot = F)
        yall <- zh$counts
        K <- length(yall)

        hist(zzz, breaks = breaks, xlab = "z-scores", main = "Empirical distribution of z-scores", freq = F)
        xfit <- seq(min(zzz), max(zzz), length = 4000)
        yfit <- dnorm(xfit / mlests[3], mean = 0, sd = 1)
        lines(xfit, yfit, col = "darkgreen", lwd = 2)
    }

    FDR <- p.adjust(pval_empirical, method = "BH")
    newList <- list("pval" = pval_empirical, "FDR" = FDR)
    return(newList)
}

#' Test function to obtain a top list of transcripts that are differentially used
#' in the contrast of interest
#'
#' @description Function to test for differential transcript usage (DTU)
#'
#' @param object A `SummarizedExperiment` instance containing a list of objects of
#'       the `StatModel` class as obtained by the `fitDTU` function of the `satuRn` package.
#'
#' @param contrasts `numeric` matrix specifying one or more contrasts of
#'        the linear model coefficients to be tested.
#'        The rownames of the matrix should be equal to the names
#'        of parameters of the model that are involved in the contrast.
#'        The column names of the matrix will be used to construct names to store
#'        the results in the rowData of the SummarizedExperiment.
#'
#' @param plot `boolean(1)` Logical, defaults to FALSE. If set to TRUE,
#'      a plot of the histogram of the empirical z-scores and the standard normal
#'      distribution will bee displayed.
#'
#' @param sort `boolean(1)` Logical, defaults to FALSE. If set to TRUE, the output
#       output of the topTable test function is sorted according to the empirical p-values.
#'
#' @examples
#' # TODO
#' @return An updated `SummarizedExperiment` that contains the `Dataframes` displaying
#'      the significance of DTU for each transcript in each contrast of interest.
#'
#' @rdname testDTU
#'
#' @author Jeroen Gilis
#'
#' @importFrom locfdr locfdr
#'
#' @export

testDTU <- function(object, contrasts, plot = FALSE, sort = FALSE) {
    models <- rowData(object)[["fitDTUModels"]] # call rowData only once

    for (i in 1:ncol(contrasts)) {
        estimates <- sapply(models, satuRn:::getEstimates, contrast = contrasts[, i])
        se <- sqrt(sapply(models, satuRn:::varContrast, contrast = contrasts[, i]))

        df <- sapply(models, satuRn:::getDfPosterior)
        df[!(sapply(df, length))] <- NA # replace numeric(0) with NA
        df <- unlist(df)

        t <- estimates / se
        pval <- pt(-abs(t), df) * 2
        regular_FDR <- p.adjust(pval, method = "BH") # regular FDR correction

        empirical <- p.adjust_empirical(pval, t, plot = plot) # empirical FDR correction
        empirical_pval <- empirical$pval
        empirical_FDR <- empirical$FDR

        result_contrast <- data.frame(
            estimates, se, df, t, pval, regular_FDR, empirical_pval,
            empirical_FDR
        )

        if (sort == TRUE) {
            result_contrast <- result_contrast[order(result_contrast$empirical_pval), ]
        }
        rowData(object)[[paste0("fitDTUResult_", colnames(contrasts)[i])]] <- result_contrast
    }
    return(object)
}

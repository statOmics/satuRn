# Compute usage estimates (log-odds) in a contrast of interest
getEstimates <- function(object, contrast) {
    coef <- getCoef(object)
    if (is.null(coef)) {
        coef <- rep(NA, times = length(contrast))
    }
    return(contrast %*% coef)
}

# Compute the variance on the usage estimates in a contrast of interest
varContrast <- function(object, contrast) {
    if (object@type != "fitError") {
        if (nrow(object@params$vcovUnsc) == length(contrast)) {
            vcovTmp <- object@params$vcovUnsc * object@varPosterior
            return(diag(t(contrast) %*% vcovTmp %*% contrast))
        }
    }
    return(NA)
}

# The code is based on the source code of the unexported loccov function
# the locFDR package from CRAN https://CRAN.R-project.org/package=locfdr.
# The code is identical to locfdr:::loccov.
locfdr_loccov <- function (N, N0, p0, d, s, x, X, f, JV, Y, i0, H, h, sw) {
    
    M = rbind(1, x - Y[1], x^2 - Y[2])
    if (sw == 2) {
        K = length(x)
        K0 = length(i0)
        toprow = c(1 - N0/N, -t(h) %*% JV/s)
        botrow = cbind(0, JV/p0)
        mat = rbind(toprow, botrow)
        M0 = M[, i0]
        dpds.dy0 = mat %*% M0/N/H[1]
        dy0.dy = matrix(0, K0, K)
        dy0.dy[, i0] = diag(1, K0)
        dpds.dy = dpds.dy0 %*% dy0.dy
        rownames(dpds.dy) = c("p", "d", "s")
        return(dpds.dy)
    }
    else {
        xstd = (x - d)/s
        U = cbind(xstd - H[2]/H[1], xstd^2 - H[3]/H[1])
        M[, -i0] = 0
        dl0plus.dy = cbind(1 - N0/N, U %*% JV/s) %*% M/N/H[1]/p0
        G <- t(X) %*% (f * X)
        dl.dy = X %*% solve(G) %*% t(X)
        dlfdr.dy = dl0plus.dy - dl.dy
        if (sw == 3) 
            return(dlfdr.dy)
        else {
            Cov.lfdr = dlfdr.dy %*% (f * t(dlfdr.dy))
            return(Cov.lfdr)
        }
    }
}

# The code is based on the source code of the unexported locmle function
# the locFDR package from CRAN https://CRAN.R-project.org/package=locfdr.
# The code is almost identical to locfdr:::locmle, with exception of some minor
# adaptations to the coding style which were required to abide by 
# the current Bioconductor coding style guides (2021-02-25)
locfdr_locmle <- function(z, xlim, Jmle = 35, d = 0, s = 1, 
                            ep = 1/1e+05, sw = 0, Cov.in) {
    
    N = length(z)
    if (missing(xlim)) {
        if (N > 5e+05) 
            b = 1
        else b = 4.3 * exp(-0.26 * log(N, 10))
        xlim = c(median(z), b * diff(quantile(z)[c(2, 4)])/(2 * qnorm(0.75)))
    }
    aorig = xlim[1] - xlim[2]
    borig = xlim[1] + xlim[2]
    z0 = z[which(z >= aorig & z <= borig)]
    N0 = length(z0)
    Y = c(mean(z0), mean(z0^2))
    that = N0/N
    for (j in seq_len(Jmle)) {
        bet = c(d/s^2, -1/(2 * s^2))
        aa = (aorig - d)/s
        bb = (borig - d)/s
        H0 = pnorm(bb) - pnorm(aa)
        fa = dnorm(aa)
        fb = dnorm(bb)
        H1 = fa - fb
        H2 = H0 + aa * fa - bb * fb
        H3 = (2 + aa^2) * fa - (2 + bb^2) * fb
        H4 = 3 * H0 + (3 * aa + aa^3) * fa - (3 * bb + bb^3) * 
            fb
        H = c(H0, H1, H2, H3, H4)
        r = d/s
        I = matrix(rep(0, 25), 5)
        for (i in seq(from=0,to=4)) I[i + 1, 0:(i + 1)] = choose(i, 0:i)
        u1 = s^(0:4)
        II = pmax(row(I) - col(I), 0)
        II = r^II
        I = u1 * (I * II)
        E = as.vector(I %*% H)/H0
        E1 = E[2]
        E2 = E[3]
        E3 = E[4]
        E4 = E[5]
        mu = c(E1, E2)
        V = matrix(c(E2 - E1^2, E3 - E1 * E2, E3 - E1 * E2, E4 - 
                         E2^2), 2)
        bett = bet + solve(V, Y - mu)/(1 + 1/j^2)
        if (bett[2] > 0) 
            bett = bet + 0.1 * solve(V, Y - mu)/(1 + 1/j^2)
        if (is.na(bett[2])) 
            break
        else if (bett[2] >= 0) 
            break
        d = -bett[1]/(2 * bett[2])
        s = 1/sqrt(-2 * bett[2])
        if (sum((bett - bet)^2)^0.5 < ep) 
            break
    }
    if (is.na(bett[2])) {
        mle = rep(NA, 6)
        Cov.lfdr = NA
        Cor = matrix(NA, 3, 3)
    }
    else if (bett[2] >= 0) {
        mle = rep(NA, 6)
        Cov.lfdr = Cov.out = NA
        Cor = matrix(NA, 3, 3)
    }
    else {
        aa = (aorig - d)/s
        bb = (borig - d)/s
        H0 = pnorm(bb) - pnorm(aa)
        p0 = that/H0
        J = s^2 * matrix(c(1, 0, 2 * d, s), 2)
        JV = J %*% solve(V)
        JVJ = JV %*% t(J)
        mat2 = cbind(0, JVJ/N0)
        mat1 = c((p0 * H0 * (1 - p0 * H0))/N, 0, 0)
        mat = rbind(mat1, mat2)
        h = c(H1/H0, (H2 - H0)/H0)
        matt = c(1/H0, -(p0/s) * t(h))
        matt = rbind(matt, cbind(0, diag(2)))
        C = matt %*% (mat %*% t(matt))
        mle = c(p0, d, s, diag(C)^0.5)
        if (sw == 1) {
            sd = mle[4:6]
            Co = C/outer(sd, sd)
            dimnames(Co) = list(c("p0", "d", "s"), c("p0", "d", "s"))
            Cor = Co[c(2, 3, 1), c(2, 3, 1)]
        }
        if (!missing(Cov.in)) {
            i0 = which(Cov.in$x > aa & Cov.in$x < bb)
            Cov.out = locfdr_loccov(N, N0, p0, d, s, Cov.in$x, Cov.in$X, 
                                    Cov.in$f, JV, Y, i0, H, h, Cov.in$sw)
        }
    }
    names(mle) = c("p0", "del0", "sig0", "sd.p0", "sd.del0", "sd.sig0")
    mle = mle[c(2, 3, 1, 5, 6, 4)]
    out = list(mle = mle)
    if (sw == 1) {
        Cor = list(Cor = Cor)
        out = c(out, Cor)
    }
    if (!missing(Cov.in)) {
        if (Cov.in$sw == 2) {
            pds. = list(pds. = Cov.out)
            out = c(out, pds.)
        }
        else if (Cov.in$sw == 3) {
            Ilfdr = list(Ilfdr = Cov.out)
            out = c(out, Ilfdr)
        }
        else {
            Cov.lfdr = list(Cov.lfdr = Cov.out)
            out = c(out, Cov.lfdr)
        }
    }
    if ((sw == 1) | !missing(Cov.in)) 
        return(out)
    else return(mle)
}

# Compute p-values under an empirical null distribution
# The code is based on the source code of 
# the locFDR package https://CRAN.R-project.org/package=locfdr
p.adjust_empirical <- function(pvalues, tvalues, plot = FALSE) {
    zvalues <- qnorm(pvalues / 2) * sign(tvalues)
    
    # to avoid locFDR (numerical) fit errors if abs(z) is extremely large
    zvalues_mid <- zvalues[abs(zvalues) < 10]
    zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)] 

    #### start code from locFDR package
    N <- length(zvalues_mid)
    b <- 4.3 * exp(-0.26 * log(N, 10))
    med <- median(zvalues_mid)
    sc <- diff(quantile(zvalues_mid)[c(2, 4)]) / (2 * qnorm(.75))
    # initial MLEs of empirical null parameters
    mlests <- locfdr_locmle(zvalues_mid, xlim = c(med, b * sc))

    lo <- min(zvalues_mid)
    up <- max(zvalues_mid)
    bre <- 120
    breaks <- seq(lo, up, length = bre)
    zzz <- pmax(pmin(zvalues_mid, up), lo)
    zh <- hist(zzz, breaks = breaks, plot = FALSE)
    x <- (breaks[-1] + breaks[-length(breaks)]) / 2
    sw <- 0
    X <- cbind(1, poly(x, df = 7))
    zh <- hist(zzz, breaks = breaks, plot = FALSE)
    y <- zh$counts
    f <- glm(y ~ poly(x, df = 7), poisson)$fit

    Cov.in <- list(x = x, X = X, f = f, sw = sw)
    ml.out <- locfdr_locmle(zvalues_mid,
        xlim = c(mlests[1], b * mlests[2]),
        d = mlests[1], s = mlests[2], Cov.in = Cov.in
    ) # updated MLEs of empirical null parameters
    mlests <- ml.out$mle
    #### end code from locFDR package

    zval_empirical <- (zvalues - mlests[1]) / mlests[2]
    # p-values under empirical null
    pval_empirical <- 2 * pnorm(-abs(zval_empirical), mean = 0, sd = 1) 

    if (plot) {
        zval_empirical <- zval_empirical[!is.na(zval_empirical)]
        lo <- min(zval_empirical)
        up <- max(zval_empirical)

        lo <- min(lo, -1 * up) # to center the figure
        up <- max(up, -1 * lo)

        bre <- 120
        breaks <- seq(lo, up, length = bre)
        zzz <- pmax(pmin(zval_empirical, up), lo)
        zh <- hist(zzz, breaks = breaks, plot = FALSE)
        yall <- zh$counts
        K <- length(yall)

        hist(zzz, breaks = breaks, xlab = "z-scores", 
            main = "Empirical distribution of z-scores", freq = FALSE)
        xfit <- seq(min(zzz), max(zzz), length = 4000)
        yfit <- dnorm(xfit / mlests[3], mean = 0, sd = 1)
        lines(xfit, yfit, col = "darkgreen", lwd = 2)
    }

    FDR <- p.adjust(pval_empirical, method = "BH")
    newList <- list("pval" = pval_empirical, "FDR" = FDR)
    return(newList)
}

#' Test function to obtain a top list of transcripts that are
#' differentially used in the contrast of interest
#'
#' @description Function to test for differential transcript usage (DTU)
#'
#' @param object A `SummarizedExperiment` instance containing a list of objects 
#'     of the `StatModel` class as obtained by the `fitDTU` function of the 
#'     `satuRn` package.
#'
#' @param contrasts `numeric` matrix specifying one or more contrasts of 
#'     the linear model coefficients to be tested. The rownames of the matrix 
#'     should be equal to the names of parameters of the model that are 
#'     involved in the contrast. The column names of the matrix will be used 
#'     to construct names to store the results in the rowData of 
#'     the SummarizedExperiment.
#'
#' @param plot `boolean(1)` Logical, defaults to FALSE. If set to TRUE,
#'     a plot of the histogram of the empirical z-scores and the standard normal
#'     distribution will be displayed.
#'
#' @param sort `boolean(1)` Logical, defaults to FALSE. If set to TRUE, 
#'     the output of the topTable test function is sorted according to 
#'     the empirical p-values.
#'
#' @examples
#' data(sumExp_example, package = "satuRn")
#' library(SummarizedExperiment)
#' sumExp <- fitDTU(
#'     object = sumExp_example,
#'     formula = ~ 0 + group,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     verbose = TRUE
#' )
#' group <- as.factor(colData(sumExp)$group)
#' design <- model.matrix(~ 0 + group)
#' colnames(design) <- levels(group)
#' L <- matrix(0, ncol = 2, nrow = ncol(design))
#' rownames(L) <- colnames(design)
#' colnames(L) <- c("Contrast1", "Contrast2")
#' L[c("VISp.L5_IT_VISp_Hsd11b1_Endou", 
#'     "ALM.L5_IT_ALM_Tnc"), 1] <- c(1, -1)
#' L[c("VISp.L5_IT_VISp_Hsd11b1_Endou", 
#'     "ALM.L5_IT_ALM_Tmem163_Dmrtb1"), 2] <- c(1, -1)
#'
#' sumExp <- testDTU(object = sumExp, contrasts = L, plot = FALSE, sort = FALSE)
#' @return An updated `SummarizedExperiment` that contains the `Dataframes` 
#'     that display the significance of DTU for each transcript 
#'     in each contrast of interest.
#'
#' @rdname testDTU
#'
#' @author Jeroen Gilis
#'
#' @importFrom locfdr locfdr
#' @importFrom stats qnorm median quantile poly glm poisson pnorm dnorm 
#' p.adjust pt
#' @importFrom graphics hist lines
#' @importFrom SummarizedExperiment colData
#'
#' @export

testDTU <- function(object, contrasts, plot = FALSE, sort = FALSE) {
    if (is.null(rowData(object)[["fitDTUModels"]])) {
        stop("fitDTUModels is empty. Did you run fitDTU first?")    
    }
    models <- rowData(object)[["fitDTUModels"]] # call rowData only once

    for (i in seq_len(ncol(contrasts))) {
        estimates <- vapply(
            X = models,
            FUN = getEstimates,
            contrast = contrasts[, i],
            FUN.VALUE = numeric(1)
        )
        se <- sqrt(vapply(
            X = models,
            FUN = varContrast,
            contrast = contrasts[, i],
            FUN.VALUE = numeric(1)
        ))
        df <- vapply(
            X = models,
            FUN = getDfPosterior,
            FUN.VALUE = numeric(1)
        )

        df[!(vapply(df, length, numeric(1)))] <- NA # replace numeric(0) with NA
        df <- unlist(df)

        t <- estimates / se
        pval <- pt(-abs(t), df) * 2
        regular_FDR <- p.adjust(pval, method = "BH") # regular FDR correction
        
        # empirical FDR correction    
        empirical <- p.adjust_empirical(pval, t, plot = plot) 
        empirical_pval <- empirical$pval
        empirical_FDR <- empirical$FDR

        result_contrast <- data.frame(
            estimates, se, df, t, pval, regular_FDR, empirical_pval,
            empirical_FDR
        )

        if (sort == TRUE) {
            result_contrast <- result_contrast[order(
                result_contrast$empirical_pval), ]
        }
        rowData(object)[[paste0("fitDTUResult_", 
                                colnames(contrasts)[i])]] <- result_contrast
    }
    return(object)
}

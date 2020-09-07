calcDispersion <- function(model,type){

  model$dispersion <- NA

  if(type != "fitError") {
    df.r <- model$df.residual
    if(df.r > 0){
      #if(any(mod$weights==0))
      #warning("observations with zero weight not used for calculating dispersion")
      model$dispersion <- sum((model$weights*model$residuals^2)[model$weights > 0])/ df.r
    }
  }
  return(model)
}

calcVcovUnscaled <- function(model,type) {

  if (!type=="fitError") {
    p1 <- 1L:model$rank
    p <- length(model$coef)
    out <- matrix(NA,p,p)
    out[p1,p1] <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    colnames(out) <- rownames(out) <- names(model$coefficients)
    model$vcovUnsc <- out
  } else {
    model$vcovUnsc <- NA
  }
  return(model)
}

getOtherCount <- function(countData, tx2gene){
  # get tx2gene in better format
  geneForEachTx <- tx2gene$gene_id[match(rownames(countData),tx2gene$isoform_id)]
  geneForEachTx <- as.character(geneForEachTx)
  stopifnot(class(geneForEachTx) %in% c("character", "factor"))
  stopifnot(length(geneForEachTx) == nrow(countData))

  forCycle <- split(1:nrow(countData), as.character(geneForEachTx))
  all <- lapply(forCycle, function(i) {
    sct <- countData[i, , drop = FALSE]
    rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[-r, , drop = FALSE])))
    rownames(rs) <- rownames(sct)
    rs
  })

  otherCount <- do.call(rbind, all)
  otherCount <- otherCount[rownames(countData), ]
  return(otherCount)
}

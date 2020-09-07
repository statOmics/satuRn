#' Test function to obtain a top list of transcripts that are differentially used
#' in the contrast of interest
#'
#' @description Test function to assesss differentially transcript usage (DTU)
#'
#' @param models A `SummarizedExperiment` containing a list of objects of the `StatModel`
#'      class as obtained by the `fitQB` function of the `qbDTU` package.
#'
#' @param contrast A `vector` or `matrix` that specifies the contrast of interest.
#'
#' @param plot `boolean(1)` Logical, defaults to FALSE. If set to TRUE,
#'      a plot of the histogram of the empirical z-scores and the standard normal
#'      distribution.
#'
#' @param sort `boolean(1)` Logical, defaults to FALSE. If set to TRUE, the output
#       output of the topTable test funcion is sorted accoring to the emprical p-values.
#'
#' @examples #TODO
#'
#' @return A `Dataframe` displaying the significance of DTU for each transcript.
#'      A `List` of `Dataframes` will be returned if multiple contrasts of interest
#'      were provided in `contrast`.
#'
#' @rdname topTable
#'
#' @author Jeroen Gilis
#'
#' @importFrom locfdr locfdr
#'
#' @export

topTable <- function(object,contrasts,plot=FALSE,sort=FALSE) {

  models <- rowData(object)[["fitQBModels"]] # call rowData only once

  for (i in 1:ncol(contrasts)) {

    estimates <- sapply(models,getEstimates,contrast=contrasts[,i])
    se <- sqrt(sapply(models,varContrast,contrast=contrasts[,i])) # uses the limma squeezed variances
    df <- sapply(models,getDfPosterior)

    t <- estimates/se

    pval <- pt(-abs(t),df)*2

    regular_FDR <- p.adjust(pval, method = "BH") # regular FDR correction
    empirical <- p.adjust_empirical(pval,t,plot=plot) # empirical FDR correction
    empirical_pval <- empirical$pval
    empirical_FDR <-  empirical$FDR

    # Not entirely valid way for lfc cut-off
    # mcp <- 0.2 # minimal change in proportion (that is considered biologically relevant)
    # t2 <- qt(empirical_pval/2,df,lower.tail = F)*sign(t)
    # estimate_2 <- t2*se
    #
    # tstat.right <- (abs(estimate_2) - mcp)/se
    # tstat.left <- (abs(estimate_2) + mcp)/se
    # empirical_pval_mcp <- pt(tstat.right, df = df, lower.tail = FALSE) +
    #   pt(tstat.left, df = df, lower.tail = FALSE)
    #
    # empirical_FDR_mcp <- p.adjust(empirical_pval_mcp, method = "BH")

    result_contrast <- data.frame(estimates,se,df,t,pval,regular_FDR,empirical_pval,
                                  empirical_FDR)

    if (sort == TRUE) {
      result_contrast <- result_contrast[order(result_contrast$empirical_pval),]
    }
    rowData(object)[[paste0("fitQBResult_",colnames(contrasts)[i])]] <- result_contrast

  }
  return(object)
}

#' A `Matrix` with transcript-level counts derived from our case study which 
#' builds on the dataset of Tasic et al. We used Salmon (V1.1.0) to quantify 
#' all L5IT cells (both for ALM and VISp tissue) from mice with a normal 
#' eye condition. From these cells, we randomly sampled 20 cells 
#' of each of the following cell types to use for this vignette;
#' L5_IT_VISp_Hsd11b1_Endou, L5_IT_ALM_Tmem163_Dmrtb1 and L5_IT_ALM_Tnc. 
#' The data has already been leniently filtered with 
#' the `filterByExpr` function of edgeR.
#'
#' @usage data(Tasic_counts_vignette)
"Tasic_counts_vignette"

#' Metadata associated with the expression matrix `Tasic_counts_vignette`. 
#' See `?Tasic_counts_vignette` for more information on the dataset.
#'
#' @usage data(Tasic_metadata_vignette)
"Tasic_metadata_vignette"

#' A `SummarizedExperiment` derived from our case study which builds on 
#' the dataset of Tasic et al. It contains the same cells as the data object 
#' used in the vignette (see `?Tasic_counts_vignette` for more information). 
#' In this SummarizedExperiment, we performed a filtering with `filterByExpr` 
#' of edgeR with more stringent than default parameter settings 
#' (min.count = 100,min.total.count = 200, large.n = 50, min.prop = 0.9) to 
#' reduced the number of retained transcripts. We used this object to create 
#' an executable example in the help files of satuRn.
#'
#' @usage data(sumExp_example)
"sumExp_example"

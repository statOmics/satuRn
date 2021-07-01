
<!-- README.md is generated from README.Rmd. Please edit that file -->

# satuRn

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/statOmics/satuRn/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/statOmics/satuRn/actions)
<!-- badges: end -->

satuRn is a highly performant and scalable method for performing
differential transcript usage analyses.

## Installation instructions

To install the current version of `satuRn` in Bioconductor, run;

``` r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager") 
}
BiocManager::install("tradeSeq")
```

To install the development version, run;

``` r
devtools::install_github("statOmics/satuRn")
```

The installation should only take a few seconds. The dependencies of the
package are listed in the DESCRIPTION file of the package.

## Issues and bug reports

Please use <https://github.com/statOmics/satuRn/issues> to submit
issues, bug reports, and comments.

## Usage

A minimal example of the different functions for `modelling`, `testing`
and `visualizing` differential transcript usage is provided.

! See the online
[vignette](https://github.com/statOmics/satuRn/blob/master/vignettes/Vignette.Rmd)
or the satuRn [website](https://statomics.github.io/satuRn/) for a more
elaborate and reproducible example.

``` r
library(satuRn)
library(SummarizedExperiment)
```

Provide a transcript expression matrix and corresponding `colData` and
`rowData`

``` r
sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = Tasic_counts_vignette),
    colData = Tasic_metadata_vignette,
    rowData = txInfo
)

# Specify design formula from colData
metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)
```

Next, we test for differential transcript usage with the function
`testDTU`. This function takes as input the `SummarizedExperiment`
object generated by `fitDTU` and a contrast matrix or vector. The latter
is used to specify the comparison(s) of interest and can either be
generated manually or automatically with the `makeContrasts` function of
the `limma` R package.

``` r
group <- as.factor(Tasic_metadata_vignette$group)
design <- model.matrix(~ 0 + group) # constructs design matrix
colnames(design) <- levels(group)
L <- limma::makeContrasts(Contrast1 = VISp.L5_IT_VISp_Hsd11b1_Endou - ALM.L5_IT_ALM_Tnc,
                          Contrast2 = VISp.L5_IT_VISp_Hsd11b1_Endou - ALM.L5_IT_ALM_Tmem163_Dmrtb1, 
                          levels = design) # constructs contrast matrix

sumExp <- satuRn::testDTU(object = sumExp, 
                          contrasts = L, 
                          plot = FALSE, 
                          sort = FALSE)
```

The test results are now saved into the `rowData` of our
`SummarizedExperiment` object under the name `fitDTUResult_` followed by
the name of the contrast of interest (i.e. the column names of the
contrast matrix). The results can be accessed as follows:

``` r
head(rowData(sumExp)[["fitDTUResult_Contrast1"]]) # first contrast
```

Finally, we may visualize the usage of select transcripts in select
groups of interest with `plotDTU`:

``` r
sumExp <- satuRn::testDTU(
    object = sumExp,
    contrasts = L,
    plot = FALSE,
    sort = FALSE
)
```

Finally, we may visualize the usage of select transcripts in select
groups of interest with `plotDTU`:

``` r
group1 <- rownames(colData(sumExp))[colData(sumExp)$group == "VISp.L5_IT_VISp_Hsd11b1_Endou"]
group2 <- rownames(colData(sumExp))[colData(sumExp)$group == "ALM.L5_IT_ALM_Tnc"]

plots <- satuRn::plotDTU(
    object = sumExp,
    contrast = "Contrast1",
    groups = list(group1, group2),
    coefficients = list(c(0, 0, 1), c(0, 1, 0)),
    summaryStat = "model",
    transcripts = c(
        "ENSMUST00000081554",
        "ENSMUST00000195963",
        "ENSMUST00000132062"
    ),
    genes = NULL,
    top.n = 6
)
```

``` r
# Example plot from our publication:
knitr::include_graphics("https://raw.githubusercontent.com/statOmics/satuRn/master/inst/figures/README-DTU_plot.png")
```

<img src="https://raw.githubusercontent.com/statOmics/satuRn/master/inst/figures/README-DTU_plot.png" width="75%" />

## Citation

Below is the citation output from using `citation('satuRn')` in R.
Please run this yourself to check for any updates on how to cite
**satuRn**.

``` r
print(citation("satuRn"), bibtex = TRUE)
```

    ## 
    ## Gilis J (2021). _Scalable Analysis of differential Transcript Usage for
    ## bulk and single-Cell RNA-sequencing applications_. doi:
    ## 10.18129/B9.bioc.satuRn (URL: https://doi.org/10.18129/B9.bioc.satuRn),
    ## https://github.com/statOmics/satuRn - R package version 1.1.0, <URL:
    ## http://www.bioconductor.org/packages/satuRn>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {Scalable Analysis of differential Transcript Usage for bulk and single-Cell RNA-sequencing applications},
    ##     author = {Jeroen Gilis},
    ##     year = {2021},
    ##     url = {http://www.bioconductor.org/packages/satuRn},
    ##     note = {https://github.com/statOmics/satuRn - R package version 1.1.0},
    ##     doi = {10.18129/B9.bioc.satuRn},
    ##   }
    ## 
    ## Gilis J, Vitting-Seerup K, Van den Berge K, Clement L (2021). "Scalable
    ## Analysis of Differential Transcript Usage for Bulk and Single-Cell
    ## RNA-sequencing Applications." _bioRxiv_. doi: 10.1101/2021.01.14.426636
    ## (URL: https://doi.org/10.1101/2021.01.14.426636), <URL:
    ## https://www.biorxiv.org/content/10.1101/2021.01.14.426636v1>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Scalable Analysis of Differential Transcript Usage for Bulk and Single-Cell RNA-sequencing Applications},
    ##     author = {Jeroen Gilis and Kristoffer Vitting-Seerup and Koen {Van den Berge} and Lieven Clement},
    ##     year = {2021},
    ##     journal = {bioRxiv},
    ##     doi = {https://doi.org/10.1101/2021.01.14.426636},
    ##     url = {https://www.biorxiv.org/content/10.1101/2021.01.14.426636v1},
    ##   }

Please note that the `satuRn` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `satuRn` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.13/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://statOmics.github.io/satuRn) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://github.com/lcolladotor/biocthis)*.

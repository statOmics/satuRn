context("StatModel construction")

test_that("StatModel construction", {
    ## Empty model
    expect_true(validObject(StatModel()))
    ## Default model type
    expect_identical(StatModel()@type, "fitError")
    ## Fully defined model
    mod <- StatModel(
        type = "glm",
        params = list(x = 3, y = 7, b = 4),
        varPosterior = c(0.1, 0.2, 0.3),
        dfPosterior = c(6, 7, 8)
    )
    expect_true(validObject(mod))
})

##############################################################

context("Statmodel accession") # 2/5 accessors tested

data("sumExp_example", package = "satuRn")

sumExp_example <- fitDTU(
    object = sumExp_example,
    formula = ~ 0 + group,
    parallel = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    verbose = TRUE
)

mod1 <- rowData(sumExp_example)[["fitDTUModels"]]$"ENSMUST00000165123"
mod_empty <- StatModel()

test_that("model parameters can be accessed", {

    # test getModel
    expect_identical(class(getModel(mod1)), "list")
    expect_equal(length(getModel(mod1)), 4)
    expect_identical(class(getModel(mod_empty)), "list")
    expect_equal(length(getModel(mod_empty)), 0)

    # test getCoef
    expect_identical(class(getCoef(mod1)), "numeric")
    expect_identical(length(getCoef(mod1)), length(levels(as.factor(colData(sumExp_example)$group))))
})

##############################################################

context("fitDTU for different input types")

# load data
data("sumExp_example", package = "satuRn")

test_that("different fitDTU assay types give same results", {

    # matrix input type
    sumExp_mat <- sumExp_example
    sumExp_mat <- fitDTU(
        object = sumExp_mat,
        formula = ~ 0 + group,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        verbose = TRUE
    )

    # DataFrame input type
    mat <- sumExp_example@assays@data$counts
    sumExp_df <- sumExp_example
    sumExp_df@assays@data$counts <- mat
    sumExp_df <- satuRn::fitDTU(
        object = sumExp_df,
        formula = ~ 0 + group,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        verbose = TRUE
    )

    # sparseMatrix input type
    sumExp_sparse <- sumExp_example
    sumExp_sparse@assays@data$counts <- Matrix(mat, sparse = TRUE)
    sumExp_sparse <- satuRn::fitDTU(
        object = sumExp_sparse,
        formula = ~ 0 + group,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        verbose = TRUE
    )

    # DelayedMatrix input type
    sumExp_delayed <- sumExp_example
    sumExp_delayed@assays@data$counts <- DelayedArray(mat)
    sumExp_delayed <- satuRn::fitDTU(
        object = sumExp_delayed,
        formula = ~ 0 + group,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        verbose = TRUE
    )

    expect_equal(rowData(sumExp_mat)[["fitDTUModels"]], rowData(sumExp_df)[["fitDTUModels"]])
    expect_equal(rowData(sumExp_mat)[["fitDTUModels"]], rowData(sumExp_sparse)[["fitDTUModels"]])
    expect_equal(rowData(sumExp_mat)[["fitDTUModels"]], rowData(sumExp_delayed)[["fitDTUModels"]])
})

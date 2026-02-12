test_that("runCATS accepts SummarizedExperiment input", {

  skip_if_not_installed("SummarizedExperiment")

  data("gene_expression", package = "CATS")
  data("human_PPIN", package = "CATS")

  expect_true(is.data.frame(gene_expression))
  expect_gt(nrow(gene_expression), 200)

  genes <- as.character(gene_expression[[1]][1:200])

  expr  <- as.matrix(gene_expression[1:200, -1, drop = FALSE])
  rownames(expr) <- genes


  sig <- gene_expression[1][1:20,]
  de  <- gene_expression[1][21:60,]

  rownames(expr) <- genes

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = expr)
  )

  res <- runCATS(
    signature = sig,
    geneExp = se,
    PPI = human_PPIN,
    DESet = de,
    pCut = 0.01
  )

  expect_true(is.list(res) || is.character(res))
  expect_gt(length(unlist(res)), 0)
})

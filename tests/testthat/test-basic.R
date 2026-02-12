test_that("package data objects load and are well-formed", {

  data("gene_expression", package = "CATS")
  data("signature", package = "CATS")
  data("DE_genes", package = "CATS")
  data("human_PPIN", package = "CATS")

  expect_true(exists("gene_expression"))
  expect_true(exists("signature"))
  expect_true(exists("DE_genes"))
  expect_true(exists("human_PPIN"))

  expect_true(
    is.matrix(gene_expression) ||
      is.data.frame(gene_expression) ||
      inherits(gene_expression, "Matrix")
  )

  expect_true(is.data.frame(signature))
  expect_true(is.data.frame(DE_genes))

  # Extract gene names robustly
  get_genes <- function(df) {
    char_cols <- which(vapply(df, is.character, logical(1)))
    expect_gt(length(char_cols), 0)
    unique(df[[char_cols[1]]])
  }

  sig <- get_genes(signature)
  de  <- get_genes(DE_genes)

  expect_gt(length(sig), 0)
  expect_gt(length(de), 0)
})

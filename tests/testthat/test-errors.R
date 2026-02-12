test_that("input validation: obvious wrong inputs error cleanly", {
  data("gene_expression", package = "CATS")
  data("signature", package = "CATS")
  data("human_PPIN", package = "CATS")

  # signature must be a data.frame (by your API)
  expect_error(getSignatureGenes(123, gene_expression))

  # geneExp must be matrix/data.frame/SummarizedExperiment
  expect_error(getSignatureGenes(signature, 123))

  # pCut must be in [0,1] (add validation in runCATS)
  expect_error(
    runCATS(signature = signature,
            geneExp = gene_expression,
            PPI = human_PPIN,
            DESet = signature,   # any non-empty placeholder; type doesn't matter for this check
            pCut = 2)
  )
})

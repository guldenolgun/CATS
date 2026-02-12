test_that("CATS pipeline runs on a tiny deterministic subset", {

  data("gene_expression", package = "CATS")
  data("human_PPIN", package = "CATS")

  expect_true(!is.null(rownames(gene_expression)))
  expect_gt(nrow(gene_expression), 60)

  sig <- gene_expression[1][1:20,]
  de  <- gene_expression[1][21:60,]

  geneExp_small <- gene_expression[1:200,]

  s0 <- getSignatureGenes(signature = data.frame(sig), geneExp = geneExp_small)
  expect_true(is.list(s0))

  s1 <- ContextAgnosticExpansion(
    expressionList = s0,
    PPI = human_PPIN,
    pCut = 0.001,
    FCCut = 0.01
  )
  expect_true(is.list(s1))

  cor_cut <- getCorr(
    geneExp = geneExp_small,
    sizeOfSampling = 20,
    quantileCut = 0.3
  )
  expect_true(is.numeric(cor_cut))
  expect_length(cor_cut, 1)

  s2 <- ContextSpecificExpansion(
    expressionList = s1,
    corCutOff = as.numeric(cor_cut),
    FCCut = 3,
    fracCut = 0.05,
    nThread = 1
  )
  expect_true(is.character(s2) || is.list(s2))

  out <- getDEgenes(geneList = s2, DESet = de)
  expect_gt(length(out), 0)
})

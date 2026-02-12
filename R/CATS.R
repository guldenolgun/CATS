#' STEP 0: Separates the signature gene expression and the rest of the gene
#' expression
#'
#' @param signature Data frame formatted signature gene list
#' @param geneExp Gene by sample expression matrix or SummarizedExperiment
#'                  that includes both gene names and gene expression profiles.
#'                  The first column of the expression matrix should be gene
#'                  names. For SummarizedExperiment, gene names must be provided
#'                   as rownames(geneExp)
#'
#' @return Returns signature and rest of the genes' expression matrix as a
#'        list
#'
#' @examples
#' step0 <- getSignatureGenes(signature, gene_expression)
#'
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment assay
#' @importFrom stats quantile
#'
#' @export
getSignatureGenes <- function(signature, geneExp) {
    if (inherits(geneExp, "SummarizedExperiment")) {
        mat <- SummarizedExperiment::assay(geneExp)
        genes <- rownames(geneExp)
        if (is.null(genes)) genes <- rownames(mat)
        if (is.null(genes)) stop("For SummarizedExperiment, gene names must be in rownames(geneExp) or rownames(assay(geneExp)).")
        geneExp <- data.frame(gene = genes, mat, check.names = FALSE)
    }else {
        geneExp <- geneExp
    }
    if (!is.data.frame(signature)) {
        stop("`signature` must be a data.frame formatted signature gene list.")
    }
    ok_geneExp <-
        is.matrix(geneExp) ||
        is.data.frame(geneExp) ||
        inherits(geneExp, "SummarizedExperiment")

    if (!ok_geneExp) {
        stop("`geneExp` must be a matrix/data.frame or a SummarizedExperiment.")
    }
    geneExp <- na.omit(geneExp)
    signature <- data.frame(signature)
    geneExp <- data.frame(geneExp)
    common <- intersect(geneExp[, 1], signature[, 1])
    signt <- geneExp[match(common, geneExp[, 1]), ]
    notsignt <- geneExp[setdiff(seq_len(dim(geneExp)[1]),
                                match(common, geneExp[, 1])), ]
    return(list(signt, notsignt))
}

#' PREPARATION: Finds the correlation cut-off for the STEP 2
#'
#' @param geneExp Gene by sample expression matrix or SummarizedExperiment
#'                  that includes both gene names and gene expression profile.
#'                  The first column of this matrix should be gene names. For
#'                  SummarizedExperiment, gene names must be provided as
#'                  rownames(geneExp)
#' @param sizeOfSampling Sampling Size. By default, it is 2000
#' @param quantileCut The cut-off that is based on the percentile of randomly
#'                    sampled gene pairs. By default, it is 0.99.
#' @param samplingGenes Binary variable whether sampling should be based on
#'                      genes or samples.This variable should be set to the
#'                      TRUE for sampling based on genes. Otherwise, it
#'                      should be FALSE
#'
#' @return The correlation cut-off threshold based on the desired quantile
#'
#' @examples
#' corr_cut <- getCorr(geneExp = gene_expression, sizeOfSampling = 100,
#'                     quantileCut = 0.95)
#'
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom SummarizedExperiment assay
#' @export
getCorr <-
    function(geneExp,
             sizeOfSampling = 2000,
             quantileCut = 0.99,
             samplingGenes = TRUE) {

        if (inherits(geneExp, "SummarizedExperiment")) {
            genes <- rownames(geneExp)
            geneExp <- data.frame(genes, assay(geneExp))
        } else {
            geneExp <- geneExp
        }
        geneExp <- na.omit(geneExp)
        expr <- geneExp
        if (is.data.frame(expr)) expr <- expr[, -1, drop = FALSE]
        corrMat <- data.frame()
        count <- 1
        if (samplingGenes) {
            inn <- sample(seq_len(dim(expr)[1]), size = sizeOfSampling,
                          replace = FALSE)
            for (i in seq_len(sizeOfSampling / 2)) {
                for (j in ((sizeOfSampling / 2) + 1):sizeOfSampling) {
                    corrMat[count, 1] <- cor(t(expr[inn[i], ]),
                                             t(expr[inn[j], ]))
                    count <- count + 1
                }
            }
        } else{
            geneExp <- geneExp[, -1]
            if (sizeOfSampling <= dim(geneExp)[2]) {
                inn <- sample(seq_len(dim(geneExp)[2]), size = sizeOfSampling,
                              replace = FALSE)
            } else{
                inn <- sample(seq_len(dim(geneExp)[2]), size = sizeOfSampling,
                              replace = TRUE)
            }

            for (i in seq_len(sizeOfSampling / 2)) {
                for (j in ((sizeOfSampling / 2) + 1):sizeOfSampling) {
                    corrMat[count, 1] <- cor((geneExp[, inn[i]]),
                                             (geneExp[, inn[j]]))
                    count <- count + 1
                }
            }
        }

        return(quantile(corrMat[, 1], na.rm = TRUE, quantileCut))
    }

#' STEP 1: Context agnostic network expansion
#'
#' @param expressionList List of gene expression for signatures and the rest
#'                       of the genes obtained from the previous step, STEP 0
#'                       with getSignatureGenes' function.
#'
#' @param PPI Data frame formatted protein-protein interaction (PPI) networks
#'
#' @param pCut Cut-off for the fraction of genes in the signature list
#'             connected with the gene g, which is not in the signature list.
#'             By default, it is 0.05
#'
#' @param FCCut Fold change (FC) cut-off for this step. By default, it is 3
#'
#' @return Returns the expanded signature and the rest of the genes' expression
#'         matrix as a list
#'
#' @examples
#' step0 <- getSignatureGenes(signature, gene_expression)
#' step1 <- ContextAgnosticExpansion(expressionList = step0, PPI = human_PPIN,
#'                                   pCut = 0.001, FCCut = 0.2)
#'
#' @export
ContextAgnosticExpansion <-
    function(expressionList,
             PPI,
             pCut = 0.05,
             FCCut = 3) {
        signt <- expressionList[[1]]
        notsignt <- expressionList[[2]]

        if (nrow(signt) == 0) stop("Signature gene set is empty after Step 0.")

        ppi_genes <- rbind(data.frame(a = unlist(PPI[, 1])),
                           data.frame(a = unlist(PPI[, 2])))
        ppi_genes <- unique(ppi_genes)

        PPI_FC <- data.frame()
        for (i in seq_len(dim(ppi_genes)[1])) {
            if (!is.element(ppi_genes[i, 1], signt[,1])) {
                gen_Inte <- c(unlist(PPI[which(
                    unlist(PPI[, 1]) %in% ppi_genes[i, 1]), 2]),
                    unlist(PPI[which(
                        unlist(PPI[, 2]) %in% ppi_genes[i, 1]), 1]))
                aa <-  length(intersect(gen_Inte, signt[,1]))
                bb <- length(intersect(gen_Inte, notsignt[,1]))
                if ((aa / dim(signt)[1]) >= pCut) {
                    PPI_FC[i, 1] <- ((aa / dim(signt)[1])) / ((bb / length(
                        intersect(ppi_genes$a, notsignt[,1])
                    )))
                }
            }
        }

        potenGenes <- ppi_genes[which(PPI_FC$V1 >= FCCut), ]
        potenGenes <- intersect(potenGenes,
                                c(expressionList[[1]][,1],
                                  expressionList[[2]][,1]))

        index <- match(intersect(potenGenes, notsignt$gene), notsignt$gene)
        signt_new <- rbind(signt, notsignt[index, ])
        notsignt_new <- notsignt[-index, ]
        return(list(signt_new, notsignt_new))
    }

#' STEP 2: Context-specific network expansion
#'
#' @param expressionList List of gene expression for signatures and the rest
#'                       of the genes obtained from the previous step, STEP 1
#'                       with 'ContextAgnosticExpansion' function.
#'
#' @param corCutOff Alpha cut-off. It defines an edge if the correlation is
#'                  bigger than the alpha. Here we term alpha as a ratio
#'                  between the fraction of genes in Sig-PPI-Expanded
#'                  connected with g and the fraction of genes genome-wide
#'                  connected with g.
#'
#' @param FCCut fold change (FC) cut-off for this step. By default, it is 3
#'
#' @param fracCut The fraction cut-off for the genes in the previously expanded
#'                gene set that are connected with gene g which is not in the
#'                expanded gene list.
#'
#' @param nThread The number thread that will be used in the analysis
#'
#' @return Returns expanded gene signatures in the given context list
#'
#' @examples
#' step0 <- getSignatureGenes(signature, gene_expression)
#' step1 <- ContextAgnosticExpansion(step0, human_PPIN, pCut = 0.01)
#' step2 <- ContextSpecificExpansion(step1)
#'
#'
#' @import WGCNA
#' @export
ContextSpecificExpansion <-
    function(expressionList,
             corCutOff = 0.3,
             FCCut = 3,
             fracCut = 0.05,
             nThread = 5) {
        signt_new <- expressionList[[1]]
        notsignt_new <- expressionList[[2]]

        FC <- data.frame()
        for (i in seq_len(dim(notsignt_new)[1])) {
            aa <- WGCNA::cor(t(notsignt_new[i, 2:dim(notsignt_new)[2]]),
                             t(signt_new[, 2:dim(signt_new)[2]]),
                             nThreads = nThread)
            leng_aa <- length(which((aa) >= corCutOff))
            if ((leng_aa / dim(signt_new)[1]) >= fracCut) {
                bb <- WGCNA::cor(t(notsignt_new[i, 2:dim(notsignt_new)[2]]),
                                 t(notsignt_new[, 2:dim(notsignt_new)[2]]),
                                 nThreads = nThread)
                leng_bb <- length(which((bb) >= corCutOff))
                FC[i, 1] <- (leng_aa / dim(signt_new)[1]) /
                    ((leng_bb / dim(notsignt_new)[1]))
            }
        }

        expand_gene <- c(signt_new[,1], notsignt_new[,1][which(FC$V1 >= FCCut)])

        return(expand_gene)
    }

#' STEP 3: Pruning
#'
#' @param geneList List of gene expression for signatures and the rest of the
#'                       genes obtained from the previous step, STEP 2 with
#'                       'ContextSpecificExpansion' function.
#'
#' @param DESet Differentially expressed gene set or list of interest in
#'              character format
#'
#' @return Returns pruned gene list
#'
#' @examples
#' step0 <- getSignatureGenes(signature, gene_expression)
#' step1 <- ContextAgnosticExpansion(step0, human_PPIN, pCut = 0.01)
#' step2 <- ContextSpecificExpansion(step1)
#' step3 <- getDEgenes(geneList = step2, DESet = DE_genes)
#'
#' @export
getDEgenes <- function(geneList, DESet) {
    if(is.data.frame(DESet))
        return(intersect(geneList, DESet[,1][[1]]))
    else
        return(intersect(geneList, DESet))
}

#' The pipeline for the CATS.
#'
#' @param signature Data frame formatted signature gene list
#' @param geneExp Gene by sample expression matrix. The first column of this
#'                matrix should be gene names.
#' @param PPI Data frame formatted protein-protein interaction (PPI) networks
#' @param pCut Cut-off for the fraction of genes in the signature list
#'             connected with the gene g, which is not in the signature list.
#'             By default, it is 0.05. This parameter is for
#'             ContextAgnosticExpansion step.
#' @param step1_FCCut Fold change (FC) cut-off for the Step 1. By default,
#'                    it is 3
#' @param corCutOff Alpha cut-off. It defines an edge if the correlation is
#'                  bigger than the alpha. Here we term alpha as a ratio
#'                  between the fraction of genes in Sig-PPI-Expanded
#'                  connected with g and the fraction of genes genome-wide
#'                  connected with g.
#' @param step2_FCCut fold change (FC) cut-off for the step 2.
#'                    By default, it is 3
#' @param fracCut The fraction cut-off for the genes in the previously expanded
#'                gene set that are connected with gene g which is not in the
#'                expanded gene list.
#' @param nThread The number thread that will be used in the analysis
#' @param DESet Differentially expressed gene set or list of interest in
#'              character format
#'
#' @return Gene set
#'
#' @examples
#' step3 <- runCATS(signature = signature,
#'                 geneExp = gene_expression,
#'                 PPI = human_PPIN,
#'                 DESet = DE_genes,
#'                 pCut = 0.01)
#'
#' @export
runCATS <- function(signature,
                    geneExp,
                    PPI,
                    DESet,
                    pCut = 0.05,
                    step1_FCCut = 3,
                    corCutOff = 0.3,
                    step2_FCCut = 3,
                    fracCut = 0.05,
                    nThread = 5) {
    if (!is.numeric(pCut) || length(pCut) != 1 || is.na(pCut) || pCut < 0 || pCut > 1) {
        stop("`pCut` must be a single numeric value in [0, 1].")
    }

    step0 <- getSignatureGenes(signature = signature, geneExp = geneExp)
    step1 <- ContextAgnosticExpansion(expressionList = step0,
                                      PPI = PPI,
                                      pCut = pCut,
                                      FCCut = step1_FCCut)
    step2 <- ContextSpecificExpansion(
        expressionList = step1,
        corCutOff = corCutOff,
        FCCut = step2_FCCut,
        fracCut = fracCut,
        nThread = nThread)
    step3 <- getDEgenes(geneList = step2, DESet = DESet)
    return(step3)
}

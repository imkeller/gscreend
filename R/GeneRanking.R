# This method is adapted from Li et al and others
# But I have to write it from scratch because there is no
# available R package (only one on github, but this would be to difficult to
# install)


# helper functions, according to definition in Li et al.
# probability needs to be transformed by beta distribution

alphaBeta <- function(p_test) {
    p_test <- sort(p_test)
    n <- length(p_test)
    return(min(pbeta(p_test, 1:n, n - (1:n) + 1)))
}

# calculate rho value
makeRhoNull <- function(n, p, nperm) {
    mclapply(1:nperm, function(x) {
        p_test <- sort.int(sample(p, n, replace = FALSE))
        alphaBeta(p_test)
    }, mc.cores = 6)
}


#' Calculate gene rank
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
calculateGeneRank <- function(object) {
    pvals <- samplepval(object)
    genes <- rbind(data.frame(gene = rowData(object@sgRNAData)$gene),
                   data.frame(gene = rowData(object@sgRNAData)$gene))

    # all pvalues lower than the threshold are set to a score of 1
    alpha_cutoff <- object@FittingOptions$alphaCutoff
    cut.pvals <- pvals <= alpha_cutoff
    score_vals <- rank(pvals) / nrow(pvals)
    score_vals[!cut.pvals ] <- 1

    # calculate rho for every count gene
    rho <- unsplit(sapply(split(score_vals, genes), alphaBeta), genes)

    guides_per_gene <- sort(unique(table(genes)))

    # store this as model parameter
    set.seed(123)
    permutations=10 * nrow(unique(genes))

    rho_nullh <- lapply(guides_per_gene,
                          makeRhoNull,
                          score_vals,
                          permutations)

    #rho_nullh <- lapply(guides_per_gene, makeRhoNull, score_vals, permutations)
    # Split by gene, make comparison with null model from makeRhoNull, and unsplit by gene
    pvalue_gene <- mclapply(split(rho, genes), function(x) {
        n_sgrnas = length(x)
        mean(rho_nullh[
            guides_per_gene == n_sgrnas][[1]] <= x[[1]])
    }, mc.cores=6)


    fdr_gene <- p.adjust(pvalue_gene, method = "fdr")

    # build new summarized experiment for the GeneData slot
    rowData <- data.frame(gene = names(pvalue_gene))
    colData <- data.frame(samplename = c("T1"),
                          timepoint = c("T1"))

    object@GeneData <- SummarizedExperiment(assays=list(
                                     pvalue = as.matrix(pvalue_gene),
                                     fdr = as.matrix(fdr_gene)),
                         rowData=rowData, colData=colData)
    object

}


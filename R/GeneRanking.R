# This method is adapted from Li et al. and others But I have
# to write it from scratch because there is no available R package
# on CRAN or Bioconductor

# helper functions, according to definition in Li et al.
# probability needs to be transformed by beta distribution

alphaBeta <- function(p_test) {
    p_test <- sort(p_test)
    n <- length(p_test)
    return(min(stats::pbeta(p_test, seq_len(n), n - seq_len(n) + 1)))
}

# calculate rho value
makeRhoNull <- function(n, p, nperm, n_cores) {
    rhonull <- BiocParallel::bplapply(seq_len(nperm), function(x) {
        p_test <- sort.int(sample(p, n, replace = FALSE))
        alphaBeta(p_test)
    })
    unlist(rhonull)
}

#' @import methods
#'
calculateGenePval <- function(pvals, genes, alpha_cutoff, n_cores) {
    cut.pvals <- pvals <= alpha_cutoff
    # ranking and scoring according to pvalues
    score_vals <- rank(pvals)/length(pvals)
    score_vals[!cut.pvals] <- 1

    # calculate rho for every count gene
    rho <- unsplit(vapply(split(score_vals, genes),
                                FUN = alphaBeta,
                                FUN.VALUE = numeric(1)),
                                genes)

    guides_per_gene <- sort(unique(table(genes)))

    # store this as model parameter
    permutations = 10 * nrow(unique(genes))

    # this does not need to be parallelized because its calling
    # a function that is already serialized
    rho_nullh <- vapply(guides_per_gene,
                        FUN = makeRhoNull,
                        p = score_vals,
                        nperm = permutations,
                        n_cores = n_cores,
                        FUN.VALUE = numeric(permutations))

    # Split by gene, make comparison with null model
    # from makeRhoNull, and unsplit by gene

    pvalue_gene <- BiocParallel::bplapply(split(rho, genes), function(x) {
        n_sgrnas = length(x)
        mean(rho_nullh[, guides_per_gene == n_sgrnas] <= x[[1]])
    })

    pvalue_gene
}


calculateGeneLFC <- function(lfcs_sgRNAs, genes) {
    # Gena LFC : mean LFC of sgRNAs
    vapply(split(lfcs_sgRNAs, genes), FUN = mean, FUN.VALUE = numeric(1))
}

#' Calculate gene rank
#'
#' @param object PoolScreenExp object
#' @param alpha_cutoff alpha cutoff for alpha-RRA (default: 0.05)
#' @param n_cores number of cores to be used (default: 1)
#'
#' @return object
#'

assignGeneData <- function(object, alpha_cutoff, n_cores) {
    # p-values for neg LFC were calculated from model
    pvals_neg <- samplepval(object)
    # p-values for pos LFC: 1 - neg.pval
    pvals_pos <- 1 - samplepval(object)

    # genes (append gene list as many times as replicates)
    n_repl <- dim(pvals_neg)[2]
    genes <- do.call("rbind", replicate(n_repl,
                        data.frame(gene = rowData(object@sgRNAData)$gene),
                        simplify = FALSE))

    # calculate pvalues
    gene_pval_neg <- calculateGenePval(pvals_neg, genes, alpha_cutoff, n_cores)
    gene_pval_pos <- calculateGenePval(pvals_pos, genes, alpha_cutoff, n_cores)

    # calculate fdrs from pvalues
    fdr_gene_neg <- stats::p.adjust(gene_pval_neg, method = "fdr")
    fdr_gene_pos <- stats::p.adjust(gene_pval_pos, method = "fdr")

    # calculate gene lfc
    lfcs_sgRNAs <- samplelfc(object)
    gene_lfc <- calculateGeneLFC(lfcs_sgRNAs, genes)

    # build new summarized experiment for the GeneData slot
    # assuming that gene order is same in neg and pos
    rowData <- data.frame(gene = names(gene_pval_neg))
    colData <- data.frame(samplename = c("T1"), timepoint = c("T1"))

    # build a summarized experiment that contains p values and fdrs
    object@GeneData <- SummarizedExperiment(
        assays = list(pvalue_neg = as.matrix(gene_pval_neg),
            fdr_neg = as.matrix(fdr_gene_neg),
            pvalue_pos = as.matrix(gene_pval_pos),
            fdr_pos = as.matrix(fdr_gene_pos),
            lfc = as.matrix(gene_lfc)),
        rowData = rowData, colData = colData)
    object

}

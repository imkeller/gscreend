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


calculateGenePval <- function(pvals, genes, alpha_cutoff) {
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
    
    # Split by gene, make comparison with null model from makeRhoNull,
    # and unsplit by gene
    pvalue_gene <- mclapply(split(rho, genes), function(x) {
        n_sgrnas = length(x)
        mean(rho_nullh[
            guides_per_gene == n_sgrnas][[1]] <= x[[1]])
    }, mc.cores=6)
    
    pvalue_gene
}


calculateGeneLFC <- function(lfcs_sgRNAs, genes) {
    # Gena LFC : mean LFC of sgRNAs
    sapply(split(lfcs_sgRNAs, genes), mean)
}

#' Calculate gene rank
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
assignGeneData <- function(object) {
    # p-values for neg LFC were calculated from model
    pvals_neg <- samplepval(object)
    # p-values for pos LFC: 1 - neg.pval
    pvals_pos <- 1 - samplepval(object)
   
    # genes (append gene list as many times as replicates)
    n_repl <- dim(pvals_neg)[2]
    genes <- do.call("rbind", 
            replicate(n_repl, 
                data.frame(gene = rowData(object@sgRNAData)$gene), 
                simplify = FALSE))

    # all pvalues lower than the threshold are set to a score of 1
    alpha_cutoff <- object@FittingOptions$alphaCutoff
    
    # calculate pvalues
    gene_pval_neg <- calculateGenePval(pvals_neg, genes, alpha_cutoff)
    gene_pval_pos <- calculateGenePval(pvals_pos, genes, alpha_cutoff) 
    
    # calculate fdrs from pvalues
    fdr_gene_neg <- p.adjust(gene_pval_neg, method = "fdr")
    fdr_gene_pos <- p.adjust(gene_pval_pos, method = "fdr")
    
    # calculate gene lfc
    lfcs_sgRNAs <- samplelfc(object)
    gene_lfc <- calculateGeneLFC(lfcs_sgRNAs, genes)

    # build new summarized experiment for the GeneData slot
    # assuming that gene order is same in neg and pos
    rowData <- data.frame(gene = names(gene_pval_neg))
    colData <- data.frame(samplename = c("T1"),
                          timepoint = c("T1"))

    # build a summarized experiment that contains p values and fdrs
    object@GeneData <- SummarizedExperiment(assays=list(
                                     pvalue_neg = as.matrix(gene_pval_neg),
                                     fdr_neg = as.matrix(fdr_gene_neg),
                                     pvalue_pos = as.matrix(gene_pval_pos),
                                     fdr_pos = as.matrix(fdr_gene_pos),
                                     lfc = as.matrix(gene_lfc)),
                         rowData=rowData, colData=colData)
    object

}


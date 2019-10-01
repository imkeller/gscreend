# Wrapper funtion to run complete analysis
#' run gscreend
#'
#' @param object PoolScreenExp object
#' @param quant1 lower quantile for least quantile of squares regression
#' (default: 0.1)
#' @param quant2 upper quantile for least quantile of squares regression
#' (default: 0.9)
#' @param alphacutoff alpha cutoff for alpha-RRA (default: 0.05)
#' @param n_cores number of cores to be used (default: 1)
#'
#' @return object
#' @export
#'
#'
#' @examples raw_counts <- read.table(
#'                         system.file('extdata', 'simulated_counts.txt',
#'                         package = 'gscreend'),
#'                         header=TRUE)
#'
#'# Create the PoolScreenExp to be analyzed
#'counts_matrix <- cbind(raw_counts$library0, raw_counts$R0_0, raw_counts$R1_0)
#'
#'rowData <- data.frame(sgRNA_id = raw_counts$sgrna_id,
#'gene = raw_counts$Gene)
#'
#'colData <- data.frame(samplename = c('library', 'R1', 'R2'),
#'timepoint = c('T0', 'T1', 'T1'))
#'
#'library(SummarizedExperiment)
#'se <- SummarizedExperiment(assays=list(counts=counts_matrix),
#'rowData=rowData, colData=colData)
#'
#'pse <- createPoolScreenExp(se)
#'
#'# Run Analysis
#'pse_an <- RunGscreend(pse)
#'

RunGscreend <- function(object,
                        quant1 = 0.1, quant2 = 0.9,
                        alphacutoff = 0.05, n_cores = 1) {
    # normalize
    norm_pse <- normalizePoolScreenExp(object)
    # calculate fold changes
    lfc_pse <- calculateLFC(norm_pse)
    # sgRNA fitting
    fit_pse <- calculateIntervalFits(
        defineFittingIntervals(lfc_pse), quant1, quant2)
    # pvalues and rank sgRNAs
    pval_pse <- calculatePValues(fit_pse)
    # rank genes
    assignGeneData(pval_pse, alphacutoff, n_cores)
}

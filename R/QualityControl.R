#' Plot replicate correlation
#'
#' @param object PoolScreenExp object
#' @param rep1 Name of replicate 1
#' @param rep2 Name of replicate 2
#'
#' @return replicate_plot
#' @export
#'
#' @examples
#' # import a PoolScreenExp object that has been generated using RunGscreend()
#' pse_an <- readRDS(
#' system.file('extdata', 'gscreend_analysed_experiment.RData',
#' package = 'gscreend'))
#' plotReplicateCorrelation(pse_an, rep1 = 'R1', rep2 = 'R2')
#'
#'@importFrom graphics plot
plotReplicateCorrelation <- function(object, rep1 = "R1", rep2 = "R2") {
    se <- sgRNAData(object)
    if(rep1 %in% se$samplename & rep2 %in% se$samplename) {
        counts1 <- assays(se[, se$samplename == rep1])$normcounts
        counts2 <- assays(se[, se$samplename == rep2])$normcounts
        plot(log2(counts1 + 1), log2(counts2 + 1))
    }
    else {
        stop(
            "Error: One or more of the replicate names you provide is ",
            "not in the list of available replicates: ",
            paste(as.character(se$samplename)))
    }
}


#' Plot model parameters from the fitting
#'
#' @param object  PoolScreenExp object
#'
#' @return plot
#' @export
#'
#' @examples # import a PoolScreenExp object that has been generated using
#' # RunGscreend()
#' pse_an <- readRDS(
#' system.file('extdata', 'gscreend_analysed_experiment.RData',
#' package = 'gscreend'))
#' plotModelParameters(pse_an)
plotModelParameters <- function(object) {
    limits <- FittingIntervals(object)
    limits <- limits[seq_len(length(limits)) - 1]
    parameters <- LFCModelParameters(object)
    graphics::par(mfrow = c(2, 2))
    plot(limits, parameters[, 1], main = "Mean")
    plot(limits, parameters[, 2], main = "Sd")
    plot(limits, parameters[, 3], main = "Xi")
}


#' Extract a results table
#'
#' @param object  PoolScreenExp object
#' @param direction Whether the table should contain information on positive or
#' negative fold changes ['negative'| 'positive']
#'
#' @return plot
#' @export
#'
#' @examples # import a PoolScreenExp object that has been generated using
#' # RunGscreend()
#' pse_an <- readRDS(
#' system.file('extdata', 'gscreend_analysed_experiment.RData',
#' package = 'gscreend'))
#' ResultsTable(pse_an, direction = 'negative')
#'
ResultsTable <- function(object, direction = "negative") {
    if (direction == "negative") {
        genedataslot <- GeneData(object)
        data.frame(Name = rownames(assays(genedataslot)$fdr_neg),
                    fdr = assays(genedataslot)$fdr_neg,
                    pval = as.numeric(assays(genedataslot)$pvalue_neg[, 1]),
                    lfc = assays(genedataslot)$lfc)
    } else if (direction == "positive") {
        genedataslot <- GeneData(object)
        data.frame(Name = rownames(assays(genedataslot)$fdr_pos),
                    fdr = assays(genedataslot)$fdr_pos,
                    pval = as.numeric(assays(genedataslot)$pvalue_pos[, 1]),
                    lfc = assays(genedataslot)$lfc)
    } else {
        stop(direction, " is not a valid argument.",
            "Select positive or negative direction for results table.")
    }
}

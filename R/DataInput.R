
#' Create PoolScreenExp Experiment
#'
#' @param data Input data object containing gRNA level data
#' (SummarizedExperiment)
#'
#' @return object  PoolScreenExp object
#' @export
#'
#' @examples raw_counts <- read.table(
#'                         system.file("extdata", "simulated_counts.txt",
#'                         package = "gscreend"),
#'                         header=TRUE)
#'
#'counts_matrix <- cbind(raw_counts$library0, raw_counts$R0_0, raw_counts$R1_0)
#'
#'rowData <- data.frame(sgRNA_id = raw_counts$sgrna_id,
#'                      gene = raw_counts$Gene)
#'
#'colData <- data.frame(samplename = c("library", "R1", "R2"),
#'                      timepoint = c("T0", "T1", "T1"))
#'
#'library(SummarizedExperiment)
#'se <- SummarizedExperiment(assays=list(counts=counts_matrix),
#'                           rowData=rowData, colData=colData)
#'
#'# create a PoolScreenExp experiment
#'pse <- createPoolScreenExp(se)
createPoolScreenExp <- function(data) {
    if(is(data, "SummarizedExperiment")) {
        message("Creating PoolScreenExp object from a SummarizedExperiment object.")
        object <- createPoolScreenExpFromSE(data)
    } else {
        message("Error: input data is not of the suported input format: SummarizedExperiment.")
    }
    object
}

#' Create PoolScreenExp Experiment from summarized experiment
#'
#' @param data SummarizedExperiment object containing gRNA level data
#'
#' @return object
#'
createPoolScreenExpFromSE <- function(data) {
    object <- new("PoolScreenExp")
    object@sgRNAData <- data
    object@FittingOptions <- list(IntervalFraction = 0.1)
    object
}

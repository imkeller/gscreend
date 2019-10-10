#' @title Class to store pooled CRISPR screening experiment
#' @description
#' The \code{poolScreenExp} class is an S4 class used to store
#' sgRNA and gene related data as well as parameters necessary
#' for statistical model.
#'
#' @import SummarizedExperiment
#'
#' @slot sgRNAData A SummarizedExperiment containing
#' the data related to sgRNAs.
#' @slot FittingIntervals A vector defining the limits of the
#' intervals used for fitting of null model.
#' @slot LFCModelParameters A vector of parameters estimated
#' when fitting the null model.
#' @slot GeneData SummarizedExperiment containing the
#' data related to genes.
#' @slot FittingOptions A named list with options for fitting:
#' IntervalFraction - fraction of sgRNAs used in every
#' fitting interval (default 0.1),
#' alphaCutoff - alpha cutoff for alpha RRA algorithm (default: 0.05).
#' @export
#'
setClass("PoolScreenExp", slots = c(
    sgRNAData = "SummarizedExperiment",
    FittingIntervals = "vector",
    LFCModelParameters = "matrix",
    GeneData = "SummarizedExperiment",
    FittingOptions = "list")
)

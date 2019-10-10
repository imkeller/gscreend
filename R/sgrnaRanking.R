#' Calculate interval limits for splitting data into subsets
#'
#' @param object  PoolScreenExp object
#'
#' @return object
#' @keywords internal
defineFittingIntervals <- function(object) {
    # split intervals based on library counts use 10% of data
    # (parameter set in the object creation function
    # createPoolScreenExpFromSE())
    quantile_lims = seq(0, 1, object@FittingOptions$IntervalFraction)
    # use unique(): in case there are many 0 the first intervals are merged
    limits <- unique(stats::quantile(refcounts(object), quantile_lims))
    object@FittingIntervals <- as.integer(limits)
    object
}


#'  Fit paramters for skew normal distribution
#'
#' @param LFC logarithmic fold changes of gRNA counts
#' @param quant1 lower quantile for least quantile of squares regression
#' (default: 0.1)
#' @param quant2 upper quantile for least quantile of squares regression
#' (default: 0.9)
#'
#' @importFrom nloptr lbfgs
#' @importFrom fGarch dsnorm
#'
#' @return fit_skewnorm
#' @keywords internal
fit_least_quantile <- function(LFC, quant1, quant2) {
    # log likelihood function of 90% most central LFC values
    ll_skewnorm <- function(x) {
        mean = x[1]
        sd = x[2]
        xi = x[3]
        # left skew normal distribution used for fit because
        # the data is skewed towards negative LFC
        r = fGarch::dsnorm(LFC, mean, sd, xi)
        quant10 <- stats::quantile(r, quant1, na.rm = TRUE)
        quant90 <- stats::quantile(r, quant2, na.rm = TRUE)
        r_quant <- r[r >= quant10 & r <= quant90]
        # abs() here is ok?
        -sum(log(abs(r_quant)))
    }
    # non-linear optimization
    fit_skewnorm <- nloptr::lbfgs(c(0, 1, 0.5),
                            ll_skewnorm,
                            lower = c(-2, 0, -2), upper = c(2, 2, 2))
    fit_skewnorm
}


#' Calculate fit parameters for every subset of data
#'
#' @param object  PoolScreenExp object
#' @param quant1 lower quantile for least quantile of squares regression
#' (default: 0.1)
#' @param quant2 upper quantile for least quantile of squares regression
#' (default: 0.9)
#'
#' @return object
#' @keywords internal
calculateIntervalFits <- function(object, quant1, quant2) {
    # one fit is done for every intervall because the skew is more important
    # for perturbaions with low initial counts
    limits <- object@FittingIntervals

    # split counts into intervals needs to be done based on reference counts !
    counts_for_fit <- as.vector(refcounts(object))
    lfc_for_fit <- as.vector(samplelfc(object))
    ncols <- length(lfc_for_fit)/length(counts_for_fit)
    refcount_mask <- rep(counts_for_fit, ncols)

    # for every count, determine lower interval limit
    fits <- matrix(nrow = (length(limits) - 1), ncol = 3)
    for (i in seq_len(length(limits) - 1)) {
        lfc_subset <- lfc_for_fit[refcount_mask > limits[i] &
                                    refcount_mask < limits[i + 1]]
        fits[i, ] <- fit_least_quantile(lfc_subset, quant1, quant2)$par
    }
    object@LFCModelParameters <- fits
    message("Fitted null distribution.")
    object
}

#' Calculate p-values
#'
#' @param object PoolScreenExp object
#'
#' @return object
#'
#' @importFrom fGarch psnorm
#' @keywords internal
calculatePValues <- function(object) {
    # p value matrix needs the same dimensions as coutn data
    dimensions <- dim(sgRNAData(object))

    # empty matrix to be filled with pvalues
    assays(sgRNAData(object))$pval <- matrix(
        nrow = dimensions[1], ncol = dimensions[2])

    # previously determined fits and corresponding lower count limits
    fits <- object@LFCModelParameters
    limits <- object@FittingIntervals

    for (i in seq_len(length(limits) - 1)) {
        # subset data according to count at T0
        mask <- refcounts(object) >= limits[i] &
                    refcounts(object) < limits[i + 1]

        # data format has to be matrix
        mask <- matrix(mask, nrow = dimensions[1], ncol = dimensions[2])

        # assign actual p values from fit
        assays(sgRNAData(object))$pval[mask] <- fGarch::psnorm(
                                    assays(sgRNAData(object))$lfc[mask],
                                    fits[i, 1], fits[i, 2], fits[i, 3])
    }

    message("Calculated p-values at gRNA level.")
    object
}

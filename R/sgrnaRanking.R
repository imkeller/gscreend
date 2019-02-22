#' Calculate intervall limits for splitting data into subsets
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
defineFittingIntervals <- function(object) {
    # split intervals based on library counts
    # use 10% of data (parameter set in the object creation function)
    quantile_lims = seq(0,1,object@FittingOptions$IntervalFraction)
    # use unique(): in case there are many 0 the first intervals are merged
    limits <- unique(quantile(refcounts(object),quantile_lims))
    object@FittingIntervals <- as.integer(limits)
    object
}


#'  Fit paramters for skew normal distribution
#'
#' @param LFC
#'
#' @return
#' @export
#'
#' @examples
fit_least_quantile <- function(LFC) {
    # log likelihood function of 90% most central LFC values
    ll_skewnorm <- function(x) {
        mean = x[1]
        sd = x[2]
        xi = x[3]
        # left skew normal distribution used for fit
        # because the data is skewed towards negative LFC
        r = dsnorm(-LFC, mean, sd, xi)
        quant90 <- quantile(r, 0.1, na.rm=TRUE)
        r_quant <- r[r >= quant90]
        # abs() here is ok?
        -sum(log(abs(r_quant)))
    }
    # non-linear optimization
    fit_skewnorm <- lbfgs(c(0,1, 1.5), ll_skewnorm,
                            lower = c(-5, -5, -5), upper = c(5,5,5))
    fit_skewnorm
}


#' Calculate fit parameters for every subset of data
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
calculateIntervalFits <- function(object) {
    # one fit is done for every intervall because the skew is more important
    # for perturbaions with low initial counts
    limits <- object@FittingIntervals

    # split counts into intervals
    counts_for_fit <- as.vector(samplecounts(object))
    lfc_for_fit <- as.vector(samplelfc(object))

    # for every count, determine lower interval limit
    fits <- rep(NA, (length(limits)-1))
    for (i in 1:(length(limits)-1)) {
        lfc_subset <- lfc_for_fit[counts_for_fit > limits[i] &
                                            counts_for_fit < limits[i+1]]
        fits[i] <- list(fit_least_quantile(lfc_subset)$par)
    }
    object@LFCModelParameters <- fits
    object
}

#' Calculate p-values
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
calculatePValues <- function(object) {
    # p value matrix needs the same dimensions as coutn data
    dimensions <- dim(se)

    # empty matrix to be filled with pvalues
    assays(object@sgRNAData)$pval <- matrix(nrow = dimensions[1],
                              ncol = dimensions[2])

    # previously determined fits and corresponding lower count limits
    fits <- object@LFCModelParameters
    limits <- object@FittingIntervals

    for (i in 1:(length(limits)-1)) {
        # subset data according to count at T0
        mask <- refcounts(object) >= limits[i] &
                refcounts(object) < limits[i+1]

        # data format has to be matrix
        mask <- matrix(mask, nrow=dimensions[1] , ncol=dimensions[2])

        # assign actual p values from fit
        assays(object@sgRNAData)$pval[mask] <-
        psnorm(assays(object@sgRNAData)$lfc[mask],
               fits[[i]][1], fits[[i]][2], fits[[i]][3])
    }

    object
}

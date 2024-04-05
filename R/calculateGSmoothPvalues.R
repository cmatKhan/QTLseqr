#' Non-parametric estimation of the null distribution of G'
#'
#' Estimate p-values for the smoothed G statistic (sometimes called G') based
#' on the non-parametric estimation method described in Magwene et al. 2011.
#' The smoothed G statistic vector is filtered to exclude outlier regions based
#' on the `outlier_vector` -- see \code{\link{createBSAResults}} for more
#' details. The mode of the filtered G' set is estimated using the \code{\link[modeest]{mlv}}
#' with options `bw = 0.5, method = 'hsm'`. Finally, the mean and variance
#' of the set are estimated using the median and mode and p-values
#' are estimated from a log normal distribution.
#'
#' @param g_smooth_vector a vector of smoothed G values (tricube weighted G statistics, sometimes called G')
#' @param outlier_vector is a boolean vector where if a value is TRUE, the corresponding
#'   value in g_smooth_vector is considered an outlier
#'
#' @importFrom modeest mlv
#' @importFrom stats median plnorm
#'
#' @export
calculateGSmoothPvalues <- function(g_smooth_vector, outlier_vector) {

  filtered_g_smooth_vector = g_smooth_vector[!outlier_vector]

  filtered_g_smooth_stats = list(
    median = stats::median(filtered_g_smooth_vector),
    mode = modeest::mlv(x = filtered_g_smooth_vector, bw = 0.5, method = "hsm")[1])

  filtered_g_smooth_stats$muE = log(filtered_g_smooth_stats$median)
  filtered_g_smooth_stats$varE = abs(filtered_g_smooth_stats$muE - log(filtered_g_smooth_stats$mode))

  1 - stats::plnorm(q = g_smooth_vector,
                    meanlog = filtered_g_smooth_stats$muE,
                    sdlog = sqrt(filtered_g_smooth_stats$varE))
}

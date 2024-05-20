#' Calculates the G statistic
#'
#' The G statisic G is defined by the equation: \deqn{G = 2*\sum_{i=1}^{4}
#' n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))}
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference
#' and alternate allele depths for each bulk, as described in the following
#' table: \tabular{rcc}{ Allele \tab High Bulk \tab Low Bulk \cr Reference \tab
#' \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr}
#' ...and \eqn{obs(n_i)} are the observed allele depths as described in the data
#' frame. Method 1 calculates the G statistic using expected values assuming
#' read depth is equal for all alleles in both bulks: \deqn{exp(n_1) = ((n_1 +
#' n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)} \deqn{exp(n_2) = ((n_2 + n_1)*(n_2
#' + n_4))/(n_1 + n_2 + n_3 + n_4)} etc...
#'
#' @param depth_list A PopulationDepthList object containing the reference and
#'  alternate allele depths for each bulk. See \code{\link{populationDepths}}
#'
#' @return A vector of G statistic values with the same length as
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1002255}{The Statistics
#'   of Bulk Segregant Analysis Using Next Generation Sequencing}
#'
#' @export
calculateGStatistic <- function(depth_list) {
  # ensure that depth_list is an S3 object with class 'PopulationDepthList'
  if (!inherits(depth_list, "PopulationDepthList")) {
    stop("Input is not a PopulationDepthList object")
  }

  population_1_alt_depth <- depth_list$population_1$alt
  population_1_ref_depth <- depth_list$population_1$ref
  population_2_alt_depth <- depth_list$population_2$alt
  population_2_ref_depth <- depth_list$population_2$ref

  population_1_total_depth <- population_1_ref_depth + population_1_alt_depth
  population_2_total_depth <- population_2_ref_depth + population_2_alt_depth
  ref_total_depth <- population_1_ref_depth + population_2_ref_depth
  alt_total_depth <- population_1_alt_depth + population_2_alt_depth
  total_depth <- ref_total_depth + alt_total_depth

  expected <- c(
    ref_total_depth * population_1_total_depth / total_depth,
    ref_total_depth * population_2_total_depth / total_depth,
    population_1_total_depth * alt_total_depth / total_depth,
    alt_total_depth * population_2_total_depth / total_depth
  )

  observed <- c(population_1_ref_depth,
                population_2_ref_depth,
                population_1_alt_depth,
                population_2_alt_depth)

  # return the G statistic, which is a vector of length equal to the
  # number of SNPs
  2 * (rowSums(observed * log(
    matrix(observed, ncol = 4) / matrix(expected, ncol = 4)
  )))
}

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
#' @param g_smoothed_vector a vector of smoothed G values (tricube weighted G statistics, sometimes called G')
#' @param outlier_vector is a boolean vector where if a value is TRUE, the corresponding
#'   value in g_smoothed_vector is considered an outlier
#'
#' @importFrom modeest mlv
#' @importFrom stats median plnorm
#'
#' @export
calculateGSmoothPvalues <- function(g_smoothed_vector, outlier_vector) {

  filtered_g_smoothed_vector = g_smoothed_vector[!outlier_vector]

  filtered_g_smoothed_stats = list(
    median = stats::median(filtered_g_smoothed_vector),
    mode = modeest::mlv(x = filtered_g_smoothed_vector, bw = 0.5, method = "hsm")[1])

  filtered_g_smoothed_stats$muE = log(filtered_g_smoothed_stats$median)
  filtered_g_smoothed_stats$varE = abs(filtered_g_smoothed_stats$muE - log(filtered_g_smoothed_stats$mode))

  1 - stats::plnorm(q = g_smoothed_vector,
                    meanlog = filtered_g_smoothed_stats$muE,
                    sdlog = sqrt(filtered_g_smoothed_stats$varE))
}

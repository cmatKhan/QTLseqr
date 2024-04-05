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

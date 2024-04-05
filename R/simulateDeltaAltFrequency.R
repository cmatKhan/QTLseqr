#' Simulates the change in Alternate Allele Frequency between population_2
#'   and population_1
#'
#' The function simulates the change in alternate allele frequency between two
#'   populations with replication
#'
#' @param population_1_alt_frequency numeric. The alternate allele frequency for
#'   population 1.
#' @param population_2_alt_frequency numeric. The alternate allele frequency for
#'   population 2.
#'
#' @inheritParams simulateDeltaAltFrequencyCI
#'
#' @return Returns a vector of length replicates of simulated
#'   \eqn{\Delta}{Delta}(Alt Allele Frequency)
#'
#' @export
simulateDeltaAltFrequency <- function(depth,
                                      population_1_alt_frequency,
                                      population_2_alt_frequency,
                                      replications = 10000,
                                      alt_frequency_filter = 0.3) {
  alt_frequencies_sample <- list(
    population_1 = rbinom(
      n = replications,
      size = depth,
      prob = population_1_alt_frequency
    ) / depth,
    population_2 = rbinom(
      n = replications,
      size = depth,
      prob = population_2_alt_frequency
    ) / depth
  )
  delta_alt_frequency <- alt_frequencies_sample$population_2 - alt_frequencies_sample$population_1

  if (!is.null(alt_frequency_filter)) {
    delta_alt_frequency <-
      delta_alt_frequency[alt_frequencies_sample$population_2 >= alt_frequency_filter |
                            alt_frequencies_sample$population_1 >= alt_frequency_filter]
  }
  delta_alt_frequency
}

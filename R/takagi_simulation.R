#' Simulate a Population's Alternate Allele Frequency
#'
#' This function calculates the alternate allele frequency for a given number of
#'   individuals in a population, based on the specified population structure.
#'   It supports two population structures: F2 and RIL. For an F2 population,
#'   alleles are sampled with probabilities reflecting Mendelian inheritance
#'   for a single heterozygous locus. For RILs (Recombinant Inbred Lines),
#'   alleles are sampled assuming an equal probability of inheritance from
#'   either parent.
#'
#' @param n non-negative integer. The number of individuals in the population.
#'
#' @inheritParams simulateDeltaSnpIndexCI
#'
#' @return Returns a numeric value representing the calculated alternate allele
#'   frequency within the population. This value is a proportion ranging from
#'   0 to 1, where 0 indicates no presence of the alternate allele and 1
#'   indicates complete presence.
#'
#' @export
simulatePopulationSnpIndex <- function(n, population_structure) {
  validatePositiveInteger(n)

  switch(population_structure,
    "F2" = mean(sample(
      x = c(0, 0.5, 1),
      size = n,
      prob = c(1, 2, 1),
      replace = TRUE
    )),
    "RIL" = mean(sample(
      x = c(0, 1),
      size = n,
      prob = c(1, 1),
      replace = TRUE
    )),
    stop("population_structure must be either 'F2' or 'RIL'")
  )
}


#' Simulates the change in Alternate Allele Frequency between population_2
#'   and population_1
#'
#' The function simulates the change in alternate allele frequency between two
#'   populations with replication
#'
#' @param population_1_snp_index numeric. The alternate allele frequency for
#'   population 1.
#' @param population_2_snp_index numeric. The alternate allele frequency for
#'   population 2.
#'
#' @inheritParams simulateDeltaSnpIndexCI
#'
#' @return Returns a vector of length replicates of simulated
#'   \eqn{\Delta}{Delta}(Alt Allele Frequency)
#'
#' @export
simulateDeltaSnpIndex <- function(depth,
                                  population_1_snp_index,
                                  population_2_snp_index,
                                  replications = 10000,
                                  snp_index_filter = NULL) {
  alt_frequencies_sample <- list(
    population_1 = rbinom(
      n = replications,
      size = depth,
      prob = population_1_snp_index
    ) / depth,
    population_2 = rbinom(
      n = replications,
      size = depth,
      prob = population_2_snp_index
    ) / depth
  )
  delta_snp_index <- alt_frequencies_sample$population_2 - alt_frequencies_sample$population_1

  if (!is.null(snp_index_filter)) {
    delta_snp_index <-
      delta_snp_index[alt_frequencies_sample$population_2 >= snp_index_filter |
        alt_frequencies_sample$population_1 >= snp_index_filter]
  }
  delta_snp_index
}

#' Simulate \eqn{\Delta}{Delta}(Alt Allele Frequency) Confidence Intervals
#'
#' This function simulates \eqn{\Delta}{Delta}(Alt Allele Frequency) confidence
#'   intervals for two populations of individuals, given a specified population
#'   structure. It calculates alternate allele frequencies within each bulk
#'   and then simulates alternate allele frequency values across a range of
#'   depths to determine confidence intervals for the
#'   \eqn{\Delta}{Delta}(Alt Allele Frequency). This process is repeated for
#'   a specified number of replications to ensure robustness of the confidence
#'   interval estimates.
#'
#' @param population_1_n A positive integer specifying the number of
#'   individuals in the first bulk.
#' @param population_2_n A positive integer specifying the number of
#'   individuals in the second bulk. Defaults to the same number as
#'   `population_1_n`.
#' @param population_structure A string indicating the population structure of
#'   the individuals. Valid options are 'F2' for an F2 population and 'RIL' for
#'   Recombinant Inbred Lines. This parameter influences how alternate allele
#'   frequencies are simulated.
#' @param depth An integer vector specifying the read depths at which to
#'   simulate alternate allele frequency calls. Defaults to a sequence from
#'   1 to 100.
#' @param replications An integer specifying the number of bootstrap
#'   replications to perform. This affects the robustness and accuracy of the
#'   simulated confidence intervals. Defaults to 10000.
#' @param snp_index_filter A numeric value specifying the minimum alt
#'   allele frequency filter threshold. Used to exclude low alt allele
#'   frequency values from the simulation. Defaults to 0.3.
#' @param ci_lower_bounds A numeric vector of probabilities with values in [0,0.5],
#'   specifying the confidence levels for which to calculate confidence
#'   intervals. For example, c(0.05, 0.025) corresponds to the 95% and 97.5%
#'   confidence intervals. Defaults to c(0.05, 0.025).
#'
#' @return A data frame with each row representing a specified depth and each
#'   column representing a \eqn{\Delta}{Delta}(Alt Allele Frequency) threshold for the
#'   corresponding confidence interval at that depth. The first column
#'   specifies the depth, and subsequent columns are named according to the
#'   specified confidence intervals.
#'
#' @importFrom stringr str_remove
#'
#' @examples
#' simulateDeltaSnpIndexCI(
#'   population_1_n = 50,
#'   population_structure = "F2",
#'   depth = 1:10,
#'   replications = 100,
#'   filter = 0.3,
#'   intervals = c(0.05, 0.025)
#' )
#' @export
simulateDeltaSnpIndexCI <- function(population_1_n,
                                    population_structure,
                                    population_2_n = population_1_n,
                                    depth_vector = 1:100,
                                    replications = 10000,
                                    snp_index_filter = 0.3,
                                    ci_lower_bounds = c(0.05, 0.025)) {
  validatePositiveInteger(population_1_n, "population_1_n")
  validatePositiveInteger(population_2_n, "population_2_n")

  # if population_structure is not one of F2 or RIL, throw an error
  if (!population_structure %in% c("F2", "RIL")) {
    stop("population_structure must be either 'F2' or 'RIL'")
  }

  # Ensure valid input for lower bounds
  if (any(ci_lower_bounds <= 0 | ci_lower_bounds >= 0.5)) {
    stop("Lower bounds must be in the range [0, 0.5]")
  }

  # ensure sorting is from largest to smallest
  ci_lower_bounds <- sort(ci_lower_bounds, decreasing = TRUE)

  intervals <- vector("numeric", length(ci_lower_bounds) * 2)
  # iterate over only even indicies. The result will be a vector where
  # odd elements are the lower bound and even elements are the upper bound
  # of the confidence intervals, eg c(0.05, 0.95, 0.025, 0.975)
  for (i in seq_along(ci_lower_bounds)) {
    intervals[2 * i - 1] <- ci_lower_bounds[i]
    intervals[2 * i] <- 1 - ci_lower_bounds[i]
  }

  # Simulate allele frequencies for each population
  snp_index_sample <- list(
    population_1 =
      replicate(
        n = replications * 10,
        simulatePopulationSnpIndex(
          n = population_1_n,
          population_structure = population_structure
        )
      ),
    population_2 =
      replicate(
        n = replications * 10,
        simulatePopulationSnpIndex(
          n = population_2_n,
          population_structure = population_structure
        )
      )
  )

  # Simulate delta alt allele frequency for each depth for each
  # population. Use the simulated sample from each population to generate
  # a sample of delta alt allele frequencies.
  CI_list <- lapply(depth_vector, function(depth) {
    alt_frequencies_sample <- list(
      population_1 = sample(snp_index_sample$population_1,
        size = replications,
        replace = TRUE
      ),
      population_2 = sample(snp_index_sample$population_2,
        size = replications,
        replace = TRUE
      )
    )

    delta_snp_index_simulation <-
      simulateDeltaSnpIndex(
        depth = depth,
        population_1_snp_index = alt_frequencies_sample$population_1,
        population_2_snp_index = alt_frequencies_sample$population_2,
        replications = replications,
        snp_index_filter = snp_index_filter
      )

    # from the sample of delta alt allele frequencies, calculate the
    # quantiles for the specified intervals
    quantile(delta_snp_index_simulation, probs = intervals, na.rm = TRUE)
  })

  # Convert list to a data frame
  CI_df <- do.call(rbind, CI_list)
  CI_df <- data.frame(CI_df)
  names(CI_df) <- paste0("CI_", stringr::str_remove(intervals, "0\\."))
  CI_df$depth <- depth_vector

  CI_df
}

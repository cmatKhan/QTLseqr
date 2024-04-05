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
#' @inheritParams simulateDeltaAltFrequencyCI
#'
#' @return Returns a numeric value representing the calculated alternate allele
#'   frequency within the population. This value is a proportion ranging from
#'   0 to 1, where 0 indicates no presence of the alternate allele and 1
#'   indicates complete presence.
#'
#' @export
simulatePopulationAltFrequency <- function(n, population_structure) {
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

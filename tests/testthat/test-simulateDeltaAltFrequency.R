test_that("simulateDeltaAltFrequency returns valid deltas", {
  deltas <- simulateDeltaAltFrequency(
    depth = 10,
    population_1_alt_frequency = 0.5,
    population_2_alt_frequency = 0.6,
    replications = 100,
    alt_frequency_filter = 0
  )

  expect_true(is.numeric(deltas))
  expect_length(deltas, 100)

  deltas_filtered <- simulateDeltaAltFrequency(
    depth = 10,
    population_1_alt_frequency = 0.5,
    population_2_alt_frequency = 0.6,
    replications = 100,
    alt_frequency_filter = 0.55
  )

  expect_true(length(deltas_filtered) < length(deltas))
})

test_that("simulateDeltaAltFrequencyCI returns correct structure", {
  result <- simulateDeltaAltFrequencyCI(
    population_1_n = 10,
    population_2_n = 10,
    population_structure = "F2",
    depth_vector = 1:10,
    replications = 100,
    alt_frequency_filter = 0.3,
    ci_lower_bounds = c(0.05, 0.025)
  )

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 5)
  expect_equal(nrow(result), 10)
})

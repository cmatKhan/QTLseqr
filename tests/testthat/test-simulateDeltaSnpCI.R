test_that("simulateDeltaAltFrequencyCI handles invalid inputs", {
  expect_error(simulateDeltaAltFrequencyCI(0, "F2", depth_vector = 1:10, replications = 100))
  expect_error(simulateDeltaAltFrequencyCI(-1, "RIL", depth_vector = 1:10, replications = 100))
})


test_that("simulateDeltaAltFrequencyCI handles invalid inputs", {
  # set the random seed
  set.seed(314)
  expect_snapshot(simulateDeltaAltFrequencyCI(population_1_n = 100,
                                           "RIL",
                                           depth_vector = 1:10,
                                           replications = 100))
})

test_that("simulatePopulationAltFrequency returns valid frequencies", {
  freq <- simulatePopulationAltFrequency(50, "F2")
  expect_true(freq >= 0 && freq <= 1)

  freq <- simulatePopulationAltFrequency(50, "RIL")
  expect_true(freq >= 0 && freq <= 1)
})

test_that("simulatePopulationAltFrequency handles invalid inputs", {
  expect_error(simulatePopulationAltFrequency(0, "F2"))
  expect_error(simulatePopulationAltFrequency(-10, "RIL"))
})

test_that("simulatePopulationSnpIndex returns valid frequencies", {
  freq <- simulatePopulationSnpIndex(50, "F2")
  expect_true(freq >= 0 && freq <= 1)

  freq <- simulatePopulationSnpIndex(50, "RIL")
  expect_true(freq >= 0 && freq <= 1)
})

test_that("simulatePopulationSnpIndex handles invalid inputs", {
  expect_error(simulatePopulationSnpIndex(0, "F2"))
  expect_error(simulatePopulationSnpIndex(-10, "RIL"))
})

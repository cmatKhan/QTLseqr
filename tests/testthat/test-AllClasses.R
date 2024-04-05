test_that("BSAExperiment constructor works with empty data", {
  # Call BSAExperiment() and expect a warning
  expect_warning(emptyResult <- BSAExperiment(), "`assays` is null")

  # Now check that the object is of the correct class
  expect_true(inherits(emptyResult, "BSAExperiment"),
    info = "Object should be a BSAExperiment instance."
  )
})


test_that("BSAResults constructor works with empty data", {
  emptyResult <- BSAResults()

  # Check class -- note that the validator is run automatically, so
  # this does test that the validator works on empty data, also
  expect_true(inherits(emptyResult, "BSAResults"),
    info = "Object should be a BSAResults instance."
  )
})

# Test that a DepthList object is successfully created with valid input
test_that("PopulationDepthList is created with valid input", {
  valid_input <- list(
    population_1 = list(
      alt = c(1, 2, 3),
      ref = c(1, 2, 3)
    ),
    population_2 = list(
      alt = c(1, 2, 3),
      ref = c(1, 2, 3)
    )
  )
  expect_s3_class(
    createPopulationDepthList(valid_input),
    "PopulationDepthList"
  )
})

# Test that function stops for non-list input
test_that("Function stops for non-list input", {
  expect_error(createPopulationDepthList("not a list"), "Input is not a list")
})

# Test that function stops if list does not contain the correct names
test_that("Function stops if list does not contain the correct names", {
  incorrect_names <- list(
    wrong_name1 = c(1, 2, 3),
    wrong_name2 = c(1, 2, 3)
  )
  expect_error(
    createPopulationDepthList(incorrect_names),
    "The first level of names must be `population_1` and `population_2`"
  )
})

# Test that function stops if list elements are not the same length
test_that("Function stops if list elements are not the same length", {
  unequal_length <- list(
    population_1 = list(
      alt = c(1, 2, 3),
      ref = c(1, 2, 3)
    ),
    population_2 = list(
      alt = c(1, 2),
      ref = c(1, 2, 3)
    )
  )
  expect_error(
    createPopulationDepthList(unequal_length),
    paste0(
      "Input depth vectors for both populations ref/alt must ",
      "all be the same length. At least 1 is not."
    )
  )
})

# Test that function stops if list elements are not numeric
test_that("Function stops if list elements are not numeric", {
  not_numeric <- list(
    population_1 = list(
      alt = c("a", "b", "c"),
      ref = c(1, 2, 3)
    ),
    population_2 = list(
      alt = c(1, 2, 3),
      ref = c(1, 2, 3)
    )
  )
  expect_error(
    createPopulationDepthList(not_numeric),
    "Input list elements are not numeric"
  )
})

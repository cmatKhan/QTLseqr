test_that("depth getters all samples works", {
  bsae = bsae_obj_fixture()

  expect_snapshot(totalDepth(bsae))
  expect_snapshot(altDepth(bsae))
})

test_that("depth getters sample select works", {
  bsae = bsae_obj_fixture()

  expect_snapshot(totalDepth(bsae, sample_set = 'SLB0021'))
  expect_snapshot(altDepth(bsae, sample_set = 'SLB0025'))
})

test_that("snpIndex works", {
  bsae = bsae_obj_fixture()

  actual = snpIndex(bsae, sample = 'SLB0021')
  expect_snapshot(actual)
})

test_that("populationDepths works", {
  bsae = bsae_obj_fixture()

  actual = populationDepths(bsae,
                             population_1_sample = 'SLB0021',
                             population_2_sample = 'SLB0025')
  expect_snapshot(actual)
})

test_that("aggregate works", {
  bsae = bsae_obj_fixture(extended = TRUE)

  # Capture the output and test for a specific message
  actual <- expect_message(
    aggregate(bsae, by = ~condition+sac_day, FUN = rowSums),
    "You will may want to re-do the `comparisons` frame and add variables to the colData"
  )

  # Expect the snapshot of the result to match the expected snapshot
  expect_snapshot(actual)
})

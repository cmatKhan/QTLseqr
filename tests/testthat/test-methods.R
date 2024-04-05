test_that("createSampleComparisonFrame works", {

  bsae = bsae_obj_fixture(add_comparisons = TRUE)

  expect_snapshot(bsae@comparisons)

})

test_that("populationDepths works", {
  bsae = bsae_obj_fixture()

  actual = populationDepths(bsae,
                             population_1_sample = 'SLB0021',
                             population_2_sample = 'SLB0025')
  expect_snapshot(actual)
})

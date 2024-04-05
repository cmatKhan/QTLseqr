test_that("calculateGStatistic works", {
  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  bulk_depths = populationDepths(bsae, bsae@comparisons[1,1], bsae@comparisons[1,2])

  g_stats = calculateGStatistic(bulk_depths)

  expect_snapshot(g_stats)
})

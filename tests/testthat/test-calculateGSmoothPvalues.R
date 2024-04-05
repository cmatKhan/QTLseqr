test_that("calculatePvals works", {
  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  depths = populationDepths(bsae, bsae@comparisons[1,1], bsae@comparisons[1,2])

  g_stats = calculateGStatistic(depths)

  smoothed_stats = smoothAlleleFrequencyMetric(rowRanges(bsae), g_stats, window_size = 5000)

  outlier_filter_vector = hampelFilter(smoothed_stats)

  pvals = calculateGSmoothPvalues(g_stats, outlier_filter_vector)

  expect_snapshot(pvals)
})

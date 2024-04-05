test_that("hampelFilter works", {
  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  bulk_depths = populationDepths(bsae, bsae@comparisons[1,1], bsae@comparisons[1,2])

  g_stats = calculateGStatistic(bulk_depths)

  smoothed_stats = smoothAlleleFrequencyMetric(rowRanges(bsae), g_stats, window_size = 5000)

  outlier_filter_vector = hampelFilter(smoothed_stats)

  expect_snapshot(outlier_filter_vector)

  outlier_filter_vector_alternate = hampelFilter_alternate(smoothed_stats)

  expect_snapshot(outlier_filter_vector_alternate)
})

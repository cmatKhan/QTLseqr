test_that("smoothAlleleFrequencyMetric works on delta SNP metric", {
  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  delta_snp_index = deltaSnpIndex(bsae, bsae@comparisons[1,1], bsae@comparisons[1,2])

  smoothed_stats = smoothAlleleFrequencyMetric(rowRanges(bsae), delta_snp_index, window_size = 5000)

  expect_snapshot(smoothed_stats)
})

test_that("smoothAlleleFrequencyMetric works on G metric", {
  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  depths = populationDepths(bsae, bsae@comparisons[1,1], bsae@comparisons[1,2])

  g_stats = calculateGStatistic(depths)

  smoothed_stats = smoothAlleleFrequencyMetric(rowRanges(bsae), g_stats, window_size = 5000)

  expect_snapshot(smoothed_stats)
})

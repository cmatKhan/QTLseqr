test_that("createBSAResults works", {
  set.seed(314)

  bsae = bsae_obj_fixture(add_comparisons=TRUE)

  expect_snapshot(createBSAResults(1L, bsae, 'g', window_size = 5000))
})

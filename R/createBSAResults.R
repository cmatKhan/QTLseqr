#' Briefly, using the natural log of Gprime a median absolute deviation (MAD)
#' is calculated. The Gprime set is trimmed to exclude outlier regions
#' (i.e. QTL) based on Hampel's rule. An alternate method for filtering out
#' QTL is proposed using absolute delta SNP indeces greater than
#' a set threshold to filter out potential QTL. An estimation of the mode of the trimmed set
#' @export
createBSAResults <- function(comparison_index,
                             bsae,
                             filter_method,
                             window_size = 1e6,
                             delta_alt_frequency_filter = 0.1,
                             hampel_threshold_multiplier = 5.2,
                             ...) {

  args_list = list(...)

  .validateCreateBSAResultsArgs(list(comparison_index = comparison_index,
                                     bsae = bsae,
                                     filter_method = filter_method,
                                     window_size = window_size,
                                     delta_alt_frequency_filter = delta_alt_frequency_filter,
                                     hampel_threshold_multiplier = hampel_threshold_multiplier))

  population_1_sample <- bsae@comparisons[comparison_index, 1]
  population_2_sample <- bsae@comparisons[comparison_index, 2]

  delta_alt_frequency <- deltaAltFrequency(bsae, population_1_sample, population_2_sample)
  population_depths <- populationDepths(bsae, population_1_sample, population_2_sample)

  # calculate alt_frequency metrics
  delta_alt_frequency_smoothed <- smoothAlleleFrequencyMetric(rowRanges(bsae),
    as.vector(delta_alt_frequency),
    window_size = window_size
  )
  # a vector where the entry is TRUE if the original value IS AN OUTLIER
  delta_alt_frequency_filter_vector <- abs(delta_alt_frequency) >= abs(delta_alt_frequency_filter)

  if("depth_vector" %in% names(args_list)){
    message(paste0("The `depth_vector` passed in the `...` ",
                   "arguments is being overwritten. ",
                   "There is no need to pass `depth_vector` explicitely"))
  }

  min_depth_vector = pmin(assays(bsae[, population_1_sample])$DP,
                          assays(bsae[, population_2_sample])$DP)

  delta_alt_frequency_ci = simulateDeltaAltFrequencyCI(
    population_1_n = bsae@metadata$population_1_n,
    population_2_n = bsae@metadata$population_2_n,
    population_structure = bsae@metadata$population_structure,
    depth_vector = min(min_depth_vector, na.rm=TRUE):max(min_depth_vector, na.rm=TRUE)
  )

  delta_alt_frequency_results_df = dplyr::tibble(
    min_depth = as.vector(min_depth_vector),
    delta_alt_frequency = as.vector(delta_alt_frequency),
    delta_alt_frequency_smoothed = as.vector(delta_alt_frequency_smoothed)
  ) %>%
    dplyr::left_join(delta_alt_frequency_ci, by = c('min_depth' = 'depth'))

  # calculate G statistics
  g_stats <- calculateGStatistic(population_depths)
  g_smoothed <- smoothAlleleFrequencyMetric(rowRanges(bsae),
    g_stats,
    window_size = window_size
  )
  g_smoothed_filter <- hampelFilter(g_smoothed, hampel_threshold_multiplier)

  # set the filter_vector according to the filter_method
  filter_vector <-
    if (filter_method == "alt_frequency") {
      delta_alt_frequency_filter_vector
    } else if (filter_method == "g") {
      g_smoothed_filter
    } else {
      stop('Invalid filter method. Must be either "alt_frequency" or "g"')
    }

  # This produces G', a smoothed version of the G statistic
  g_smoothed_pvalues <- calculateGSmoothPvalues(g_smoothed, filter_vector)

  # calculate the pvalues and adj pvalues from the smooothed stat
  bh_adj_pvalues <- p.adjust(g_smoothed_pvalues, method = "BH")

  # if there are enough tests, produce the qvalue also
  tryCatch(
    {
      qvalues <- qvalue::qvalue(g_smoothed_pvalues)$qvalues
    },
    error = function(e) {
      qvalues <<- rep(NA, length(g_smoothed_pvalues))
    }
  )

  # collect the parameters that impact analysis
  analysis_params <- list(
    filter_method = filter_method,
    window_size = window_size,
    comparison = paste(population_1_sample,
      population_2_sample,
      sep = "-vs-"
    )
  )
  if (filter_method == "alt_frequency") {
    analysis_params$delta_alt_frequency_filter <- delta_alt_frequency_filter
  } else if (filter_method == "g") {
    analysis_params$hampel_threshold_multiplier <- hampel_threshold_multiplier
  }

  g_results_df = data.frame(
    g = g_stats,
    g_smooth = g_smoothed,
    filter = filter_vector,
    pvalue = g_smoothed_pvalues,
    bh_adj_pvalue = bh_adj_pvalues,
    qvalue = qvalues
  )

  res_df = dplyr::bind_cols(delta_alt_frequency_results_df, g_results_df) %>%
    as.data.frame()

  # Create and return BSAResults object
  new("BSAResults",
    res_df,
    rowRanges = rowRanges(bsae),
    analysisParams = analysis_params
  )
}

#' Check the input values of createBSAResults
#'
#' @keywords internal
.validateCreateBSAResultsArgs <- function(args_list) {

  msg <- NULL

  # bsae must be a BSAExperiment object
  if (!inherits(args_list$bsae, "BSAExperiment")) {
    msg <- c(msg, "The assays must contain at minimum 'DP' and 'AD'" )
  }
  # the comparisons slot must be a dataframe with nrow > 0
  if (!nrow(args_list$bsae@comparisons) > 0) {
    msg <- c(msg, "There are no comparisons in the BSAExperiment object")
  }
  # validate comparison_index
  if (!is.integer(args_list$comparison_index)) {
    msg <- c(msg, "The comparison index is not an integer")
  }
  if (args_list$comparison_index < 1) {
    msg <- c(msg, "The comparison index is less than 1")
  }
  if (args_list$comparison_index > nrow(args_list$bsae@comparisons)) {
    stop("The comparison index is greater than the number of comparisons")
  }
  if (!args_list$filter_method %in% c("alt_frequency", "g")) {
    msg <- c(msg, "The filter method must be either 'alt_frequency' or 'g'")
  }
  if (abs(args_list$delta_alt_frequency_filter) >= 0.5) {
    msg <- c(msg, "`delta_alt_frequency_filter` should be less than 0.5")
  }

  if(!is.null(msg)) {
    stop(msg)
  }
}

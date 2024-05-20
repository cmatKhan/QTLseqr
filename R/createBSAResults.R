#' Calculate marker frequency statistics for a objectxperiment object
#'
#' @param object A BSAExperiment object
#' @param comparison_index a single integer value indicating the row of
#'   comparisons(object) to utilize for the statistical comparison, or
#'   a vector of integers indicating multiple comparisons to calculate
#' @param statistics A character vector indicating which statistics to calculate.
#'   Currently 'delta_snp_index' and 'g' are valid choices
#'
#' @inheritDotParams smoothAlleleFrequencyMetric window_size
#' @inheritDotParams g_statistics_results outlier_filter_method delta_snp_index_threshold hampel_filter_threshold
#'
#' @importFrom foreach foreach %do%
#'
#' @export
createBSAResults <- function(object,
                             comparison_index = NULL,
                             statistics = c("delta_snp_index", "g"),
                             ...) {
  # verify that `object` is a BSAExperiment object
  if (!inherits(object, "BSAExperiment")) {
    stop("The `object` argument object must be a BSAExperiment object")
  }

  # verify that the `comparisons` slot of the BSAExperiment object is not
  # NULL
  if (is.null(comparisons(object))) {
    stop("The `comparisons` slot of the BSAExperiment object cannot be NULL")
  }

  if (is.null(comparison_index)) {
    comparison_index <- 1:nrow(comparisons(object))
  } else if (!all(comparison_index %in% 1:nrow(comparisons(object)))) {
    stop(
      "The comparison index must be an integer between 1 and ",
      "the maximum number of rows in `comparisons(object)`"
    )
  }

  res_list <- foreach::foreach(
    i = comparison_index
  ) %do% {
    comparison_samples <- comparisons(object)[i, ]

    res_list <- list()

    if ("delta_snp_index" %in% statistics) {
      tryCatch(
        {
          res_list$delta_snp_index <- delta_snp_index_results(
            object,
            comparison_samples$population_1,
            comparison_samples$population_2,
            ...
          )
        },
        error = function(e) {
          message(
            "Failed to calculate delta_snp_index statistics for samples: ",
            e$message
          )
        }
      )
    }

    if ("g" %in% statistics) {
      tryCatch(
        {
          res_list$g <- g_statistics_results(
            object,
            comparison_samples$population_1,
            comparison_samples$population_2,
            ...
          )
        },
        error = function(e) {
          message(
            "Failed to calculate g statistics for samples: ",
            e$message
          )
        }
      )
    }
    names(res_list) = NULL
    do.call('cbind', res_list)
  }
  comparison_names_list <- paste(comparisons(object)[[comparison_index,'population_2']],
                                 comparisons(object)[[comparison_index,'population_1']],
    sep = "-vs-"
  )

  names(res_list) <- comparison_names_list

  res_list <- S4Vectors::SimpleList(res_list)

  # Create and return BSAResults object
  BSAResults(
    assays = res_list,
    rowRanges = rowRanges(object)
  )
}

#' Create a single comparison DataFrame
#'
#' This is an internal function which creates a DataFrame of QTLseqr results.
#'   createBSAResult() should be used by general users as it returns the
#'   BSAResults object which includes the rowRanges and colData. However, it
#'   is exposed in order to provide documentation on the arguments which may
#'   be passed to createBSAResults() in order to adjust how smoothing and
#'   filtering is conducted, and what results are returned.
#'
#' @param object A objectxperiment object
#' @param population_1_sample One of the samples from object$sample which
#'   will act as the "low bulk" or base sample for delta snp index calculation.
#' @param population_2_sample One of the samples from object$sample which
#'  will act as the "high bulk" or comparison sample for delta snp index
#'  calculation.
#' @inheritParams smoothAlleleFrequencyMetric
#'
#' @importFrom dplyr tibble left_join
#' @importFrom S4Vectors DataFrame
#' @keywords internal
delta_snp_index_results <- function(object,
                                    population_1_sample,
                                    population_2_sample,
                                    ...) {
  # Calculate delta snp index
  delta_snp_index <- as.vector(deltaSnpIndex(
    object,
    population_1_sample,
    population_2_sample
  ))
  # calculate snp_index metrics
  delta_snp_index_smoothed <- smoothAlleleFrequencyMetric(
    rowRanges(object),
    delta_snp_index,
    ...
  )



  min_depth_vector <- as.vector(pmin(
    totalDepth(object, population_1_sample),
    totalDepth(object, population_2_sample)
  ))

  ci_depth_vector <- min(min_depth_vector, na.rm = TRUE):max(min_depth_vector, na.rm = TRUE)

  delta_snp_index_ci <- simulateDeltaSnpIndexCI(
    population_1_n = object@metadata$population_1_n,
    population_2_n = object@metadata$population_2_n,
    population_structure = object@metadata$population_structure,
    depth_vector = ci_depth_vector
  )

  # return the result as a DataFrame
  dplyr::tibble(
    min_depth = as.vector(min_depth_vector),
    delta_snp_index = as.vector(delta_snp_index),
    delta_snp_index_smoothed = as.vector(delta_snp_index_smoothed)
  ) %>%
    dplyr::left_join(delta_snp_index_ci, by = c("min_depth" = "depth")) %>%
    S4Vectors::DataFrame()
}

#' Calculate g statistics and return a DataFrame of results
#'
#' This is an internal function only -- not exposed in the API. However, it is
#'   accessible indirectly through createBSAResults(). This function calculates
#'   the G statistic for each SNP in the objectxperiment object and returns
#'   a DataFrame of the results.
#'
#' @param object A objectxperiment object
#' @param population_1_sample One of the samples from object$sample which
#'   will act as the "low bulk" or base sample for delta snp index calculation.
#' @param population_2_sample One of the samples from object$sample which
#'   will act as the "high bulk" or comparison sample for delta snp index
#'
#' @return A DataFrame of G statistic results
#'
#' @importFrom qvalue qvalue
#' @importFrom S4Vectors DataFrame
#' @keywords internal
g_statistics_results <- function(object,
                                 population_1_sample,
                                 population_2_sample,
                                 outlier_filter_method = "delta_snp_index",
                                 delta_snp_index_threshold = 0.1,
                                 hampel_filter_threshold = 5.2,
                                 ...) {
  if (!outlier_filter_method %in% c("delta_snp_index", "hampel")) {
    stop("The filter method must be either 'delta_snp_index' or 'hampel'")
  }

  if (abs(delta_snp_index_threshold) >= 0.5) {
    stop("`delta_snp_index_filter` should be less than 0.5")
  }

  # calculate G statistics
  g_stats <- as.vector(calculateGStatistic(populationDepths(
    object,
    population_1_sample,
    population_2_sample
  )))

  # set the filter_vector according to the filter_method
  metric_mask <-
    if (outlier_filter_method == "delta_snp_index") {
      delta_snp_index <- deltaSnpIndex(
        object,
        population_1_sample,
        population_2_sample
      )
      as.vector(abs(delta_snp_index) >= abs(delta_snp_index_threshold) | is.na(delta_snp_index))
    } else if (outlier_filter_method == "hampel") {
      stop("hampel filter not yet implemented.")
      # implement both a g_stat and delta_snp hampel filter in this case?
    } else {
      stop('Invalid filter method. Must be either "snp_index" or "hampel"')
    }

  g_smoothed <- smoothAlleleFrequencyMetric(
    rowRanges(object),
    g_stats,
    ...
  )

  # This produces G', a smoothed version of the G statistic
  g_smoothed_pvalues <- calculateGSmoothPvalues(g_smoothed, metric_mask)

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

  S4Vectors::DataFrame(
    g = g_stats,
    g_smoothed = g_smoothed,
    filter = metric_mask,
    pvalue = g_smoothed_pvalues,
    bh_adj_pvalue = bh_adj_pvalues,
    qvalue = qvalues
  )
}

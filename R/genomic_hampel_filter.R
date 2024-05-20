#' Complementary Error Function
#'
#' @note copied from PRACMA, which is distributed under a GNU v3+
#'
#' @author Copyright (c) 2015 Hans W Borchers
#'
#' @keywords internal
.erfinv <- function(y) {
  y[abs(y) > 1] <- NA
  sqrt(qchisq(abs(y), 1) / 2) * sign(y)
}

# Define kappa for MAD calculation
# See the top of this file for the .erfinv function details
.KAPPA <- 1 / (sqrt(2) * .erfinv(0.5)) # Approx. 1.482602

#' Hampel Filter for Genomic Data in GRanges Object
#'
#' Applies the Hampel filter to a metric in a GRanges object to identify
#'   outliers. The filter uses genomic windows to compute the median and
#'   the median absolute deviation (MAD) for each sample, then identifies
#'   outliers. See the MATLAB hampel filter documentation for more details:
#'   https://www.mathworks.com/help/dsp/ref/hampelfilter.html
#'
#' @param gr A GRanges object containing the genomic data.
#' @param metric_colname The name of the metric stored in the metadata columns (mcols) of the GRanges object.
#' @param hampel_threshold_multiplier A numeric value specifying the number of
#' median absolute deviations from the median to use as a threshold for filtering
#' outlier regions. Default is 5.2.
#' @param window_size The size of the genomic window used for computing the median and MAD, in base pairs (e.g., 1e6 for 1 Mb).
#'
#' @return A boolean vector of the same length as the input GRanges object where TRUE
#' indicates that the corresponding metric value is an outlier.
#'
#' @importFrom GenomicRanges mcols tileGenome
#' @importFrom IRanges subsetByOverlaps
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#'
#' @examples
#' gr <- GRanges(
#'   seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1"),
#'   ranges = IRanges(start = c(1, 2e5, 4e5, 6e5, 8e5), width = 1),
#'   metric = c(1, 2, 100, 3, 2)
#' )
#' outliers <- hampelFilterGRanges(gr, "metric", window_size = 1e6)
#' print(outliers)
#'
#' @export
hampelFilterGRanges <- function(gr,
                                metric_colname = "metric",
                                window_size = 2e6,
                                hampel_threshold_multiplier = 5.2,
                                ncpus = 1,
                                parallel_backend = NULL) {
  # Check if metric is available
  if (!metric_colname %in% colnames(GenomicRanges::mcols(gr))) {
    stop(
      "The specified metric",
      metric_colname,
      "is not found in the metadata columns of the GRanges object."
    )
  }

  # Register the parallel backend if ncpus > 1 or a custom backend is provided
  if (ncpus > 1 && is.null(parallel_backend)) {
    message(
      "Registering ", ncpus, " cpus for parallel processing. ",
      "The parallel backend will be stopped on exit of this function."
    )
    doParallel::registerDoParallel(cores = ncpus)
    on.exit(doParallel::stopImplicitCluster())
  } else if (!is.null(parallel_backend)) {
    message(
      "If the parallel backend is passed, `ncpus` is ignored. ",
      "Make sure your backend has the appropriate number of cpus, etc. ",
      "Additionally, the backend will not be stopped on exit of this ",
      "function."
    )
    doParallel::registerDoParallel(parallel_backend)
  }

  # Tile the genome into windows
  tiles <- GenomicRanges::tileGenome(seqinfo(gr), tilewidth = window_size, cut.last.tile.in.chrom = TRUE)

  gr$window_sets <- factor(S4Vectors::subjectHits(GenomicRanges::findOverlaps(gr, tiles)))

  window_list <- unique(gr$window_sets)

  foreach::foreach(
    window_set = levels(gr$window_sets),
    .combine = "c",
    .packages = c("GenomicRanges", "IRanges")
  ) %dopar% {
    current_gr_subset <- gr[gr$window_set == window_set]

    current_gr_subset_metric_vector <- GenomicRanges::mcols(current_gr_subset)[[metric_colname]]

    # if the metric values in the current subset are entirely NA, then return
    # TRUE, signifying that they are "outliers". May want to reconsider
    # this and set the value to NA instead of T/F?
    if (all(is.na(current_gr_subset_metric_vector))) {
      !logical(length(current_gr_subset))
    } else {
      # subset the gr object
      current_gr_subset_metric_vector_median <- median(current_gr_subset_metric_vector, na.rm = TRUE)

      # Compute the median and MAD
      window_median <- median(current_gr_subset_metric_vector, na.rm = TRUE)
      window_mad <- median(abs(current_gr_subset_metric_vector
                               - current_gr_subset_metric_vector_median),
                           na.rm = TRUE)

      # Define the threshold
      threshold <- hampel_threshold_multiplier * window_mad * .KAPPA

      # Identify outliers
      abs(current_gr_subset_metric_vector - window_median) > threshold
    }
  }
}

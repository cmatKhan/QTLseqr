#' Calculate tricube weighted statistics for each allele
#'
#' Uses local regression (wrapper for \code{\link[locfit]{locfit}}) to predict a
#' tricube smoothed version of the statistic supplied for each allele. This works as a
#' weighted average across neighboring alleles that accounts for Linkage
#' disequilibrium (LD) while minizing noise attributed to allele calling errors.
#' Values for neighboring alleles within the window are weighted by physical
#' distance from the focal allele.
#'
#' @param allele_row_ranges A GRangesList object where each record represents
#'   The genomic coordinate of a allele
#' @param stat A vector of statistics calculated from the alleles. This might be
#'   the G statistic, see \code{\link{calculateGStatistic}} or
#'   the alt_frequency from \code{\link{populationDepths}}, for example.
#'   The order of this vector must correspond with the order of the
#'   allele_row_ranges object.
#' @param window_size the window size (in base pairs) bracketing each allele for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed.
#' @param ... Other arguments passed to \code{\link[locfit]{locfit}} and
#'   subsequently to \code{\link[locfit]{locfit.raw}}() (or the lfproc). Usefull
#'   in cases where you get "out of vertex space warnings"; Set the maxk higher
#'   than the default 100. See \code{\link[locfit]{locfit.raw}}().
#'
#' @return Returns a vector of the weighted statistic caluculted with a tricube
#'   smoothing kernel
#'
#' @seealso \code{\link{getG}} for G statistic calculation
#' @seealso \code{\link[locfit]{locfit}} for local regression
#'
#' @export
smoothAlleleFrequencyMetric <- function(allele_row_ranges, stat, window_size = 2e6, ...) {
  # Ensure allele_row_ranges is a GRangesList object
  if (!inherits(allele_row_ranges, "GRanges")) {
    stop("allele_row_ranges must be a GRangesList object")
  }
  # Ensure stat length matches allele_row_ranges length
  if (length(stat) != length(allele_row_ranges)) {
    stop("Length of allele_stat must match length of allele_row_ranges")
  }

  # Add allele_stat as a metadata column to allele_row_ranges
  mcols(allele_row_ranges)$allele_stat <- stat

  # Initialize a vector to hold the smoothed statistics
  smoothed_stats <- numeric(length(stat))

  # Iterate over each tile and apply tricubeStat
  for (chrom in seqnames(allele_row_ranges)@values) {
    # Subset allele_row_ranges for the current chromosome
    alleles_in_chrom <- allele_row_ranges[seqnames(allele_row_ranges) == chrom]

    if (length(alleles_in_chrom) > 0) {
      # Extract positions and allele_stat values for alleles in this tile
      chrom_allele_positions <- start(alleles_in_chrom)
      chrom_allele_stat <- mcols(alleles_in_chrom)$allele_stat

      chrom_mean <- locfit::lp(chrom_allele_positions, h = window_size, deg = 0)
      smoothing_model <- locfit::locfit(chrom_allele_stat ~ chrom_mean, ...)
      # Apply tricubeStat
      smoothed_values <- stats::predict(smoothing_model, chrom_allele_positions)

      # Place smoothed values into the corresponding positions in smoothed_stats
      index_in_allele_row_ranges <- match(chrom_allele_positions, start(allele_row_ranges))
      smoothed_stats[index_in_allele_row_ranges] <- smoothed_values
    }
  }

  # return smoothed_stats
  smoothed_stats
}

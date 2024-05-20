# put parameterized max depth filter at some amount times the avg coverage
# for that sample

#' Extract alternative frequency data frame
#'
#' @param bsae_sample_subset A SummarizedExperiment object subset for a specific sample.
#' @inheritParams process_chromosome
#' @param verbose if TRUE, a message reporting the current sample will be
#'   displayed If 2, both sample and chromosome messages will be displayed.
#' @param ncpus Number of cores to use for parallel processing.
#' @inheritDotParams changepoint.np::cpt.np
#'
#' @return A tibble containing the processed alternative frequency data.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom SummarizedExperiment seqnames
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#'
#' @export
extract_alt_freq_dataframe <- function(bsae_sample_subset,
                                       depth_threshold = 10,
                                       ref_freq_thres = 0.2,
                                       hetero_freq_thres = c(0.4, 0.6),
                                       alt_freq_thres = 0.9,
                                       neighborhood_size = 20,
                                       verbose = FALSE,
                                       ncpus = 1,
                                       ...) {
  if (verbose) {
    message("Starting analysis on sample: ", bsae_sample_subset$sample[1])
  }

  doParallel::registerDoParallel(cores = ncpus)

  results <- foreach::foreach(
    chr = as.character(SummarizedExperiment::seqnames(bsae_sample_subset)@values),
    .combine = 'c',  # Combine results using concatenation
    .init = list(),  # Initialize with an empty list
    .multicombine = TRUE,  # Efficient combination of large results
    .packages = c("SummarizedExperiment", "dplyr")
  ) %dopar% {
    if (verbose) {
      message("Processing chromosome: ", chr)
    }
    chr_sample_subset <- bsae_sample_subset[SummarizedExperiment::seqnames(bsae_sample_subset) == chr, ]

    # step 1: calculate the allele frequency given certain thresholds and
    # smoothing parameters
    allele_df <- process_chromosome(
      chr_sample_subset, depth_threshold, ref_freq_thres,
      hetero_freq_thres, alt_freq_thres, neighborhood_size, verbose
    )

    # step 2: using the smoothed, denoised allele frequency, call breakpoints using
    # changepoint.np
    breakpoint_df <- tryCatch(
      {
        find_breakpoints(allele_df, ...)
      },
      error = function(e) {
        message("Error finding breakpoints for chromosome ", chr, ": ", e$message)
        data.frame()
      }
    )

    list(list(allele_df = allele_df, breakpoint_df = breakpoint_df))
  }

  doParallel::stopImplicitCluster()

  # bind all chromosome tables together for a given sample
  allele_dfs <- purrr::map(results, ~.$allele_df) %>% dplyr::bind_rows()
  breakpoint_dfs <- purrr::map(results, ~.$breakpoint_df) %>% dplyr::bind_rows()

  list(allele_df = allele_dfs, breakpoint_df = breakpoint_dfs)
}


#' Process Chromosome Data
#'
#' A helper function used by `extract_alt_freq_dataframe` to process data for'
#'   each chromosome separately. This function is designed to be used
#'   internally and is not exported.
#'
#' @param chr Character string specifying the chromosome name.
#' @param bsae_sample_subset A `SummarizedExperiment` object subset for a
#'   specific sample.
#' @param depth_threshold Integer specifying the minimum depth required for a
#'   variant to be considered. Default is 10.
#' @param ref_freq_thres Numeric specifying the the maximum value for which
#'   a given variant will be 'rounded' to REF (0) for the smoothed allele
#'   frequency column
#' @param hetero_freq_thres Numeric vector specifying the range of values for
#'    which a given variant will be 'rounded' to HET (0.5) for the smoothed
#'    allele frequency column
#' @param alt_freq_thres Numeric specifying the minimum value for which a given
#'   variant will be 'rounded' to ALT (1) for the smoothed allele frequency
#'   column
#'
#' @param neighborhood_size Integer specifying the number of data points around
#'   a given variant that will be considered when determining if the variant
#'   is 'noisy', or not matching its immediate neighbors. Default is 20.
#' @param verbose Logical indicating whether to print processing messages.
#'
#' @importFrom SummarizedExperiment assays seqnames start rowRanges
#' @importFrom dplyr filter mutate case_when mutate
#'
#' @return A data frame with processed allele frequency data for the specified
#'   chromosome.
process_chromosome <- function(chr_sample_subset,
                               depth_threshold,
                               ref_freq_thres,
                               hetero_freq_thres,
                               alt_freq_thres,
                               neighborhood_size,
                               verbose) {

  local_sample <- unique(chr_sample_subset$sample)
  local_seqname <- unique(as.character(SummarizedExperiment::rowRanges(chr_sample_subset)@seqnames))

  if (length(sample) > 1) {
    stop("The input data frame must have only one sample")
  }
  if (length(local_seqname) > 1) {
    stop("The input data frame must have only one seqnames")
  }

  if (verbose) {
    message("Processing chromosome: ", local_seqname)
  }

  alt_freq <- as.vector(SummarizedExperiment::assays(chr_sample_subset)$AD /
    SummarizedExperiment::assays(chr_sample_subset)$DP)
  valid_indices <- which(!is.na(alt_freq))

  if (length(valid_indices) == 0) {
    return(data.frame())
  } # Handle empty datasets early

  # step 1: extract only NA positions on a given chr for a given sample
  na_free_alt_freq <- data.frame(
    sample = local_sample,
    seqnames = local_seqname,
    pos = SummarizedExperiment::start(chr_sample_subset)[valid_indices],
    alt_freq = alt_freq[valid_indices],
    orig_index = valid_indices,
    depth = SummarizedExperiment::assays(chr_sample_subset)$DP[valid_indices]
  ) %>%
    dplyr::filter(depth >= depth_threshold)

  if (nrow(na_free_alt_freq) == 0) {
    return(data.frame())
  } # Check for no data after depth filtering

  # step 2: 'round' or smooth the allele frequency data to 0, 0.5 or 1 depending
  # on the thresholds
  rounded_data <- na_free_alt_freq %>%
    dplyr::mutate(smoothed_alt_freq = case_when(
      alt_freq < ref_freq_thres ~ 0,
      alt_freq >= hetero_freq_thres[1] & alt_freq <= hetero_freq_thres[2] ~ 0.5,
      alt_freq > alt_freq_thres ~ 1,
      TRUE ~ NA_real_
    )) %>%
    dplyr::filter(!is.na(smoothed_alt_freq))

  if (nrow(rounded_data) == 0) {
    return(data.frame())
  } # Check for no data after rounding

  # step 3: take the mode in the neighborhood around a given variant. If the
  # variant doesn't match the mode of a window of length neighborhood_size on
  # its left AND its right, it is labelled as noisy
  rounded_data$noisy_flag <- sapply(1:nrow(rounded_data), function(i) {
    is_noisy(rounded_data$smoothed_alt_freq, i, neighborhood_size)
  })

  rounded_data
}

#' Find Breakpoints in Allele Frequency Data
#'
#' Identifies breakpoints in processed allele frequency data, typically
#'   to detect genomic regions of interest such as copy number variations or
#'   significant shifts in allele frequency.
#'
#' @param df A data frame containing smoothed allele frequencies and other
#'   relevant data.
#' @inheritParams changepoint.np::cpt.np
#'
#' @importFrom changepoint.np cpt.np
#' @importFrom dplyr filter pull mutate select
#'
#' @return A data frame indicating the positions of identified breakpoints and their characteristics.
#' @export
find_breakpoints <- function(df, ...) {
  local_sample <- unique(df$sample)
  local_seqnames <- unique(df$seqnames)

  if (length(local_sample) > 1 | length(local_seqnames) > 1) {
    stop("The input data frame must have only one sample and one seqnames")
  }

  denoised_df = df %>%
    dplyr::filter(!noisy_flag)

  bp = tryCatch({
    denoised_df %>%
      dplyr::pull(smoothed_alt_freq) %>%
      changepoint.np::cpt.np(...)
    }, error = function(e){
      message('Error finding breakpoints: ', e$message)
    })

  # note, when extracting the cpts, the last index of the original dataframe
  # is always included. -1 removes it
  tryCatch({
    bp_indicies <- sort(c(bp@cpts[-length(bp@cpts)], pmin(bp@cpts[-length(bp@cpts)] + 1, max(df$pos))))
  }, error = function(e){
    message('here!')
  })
  denoised_df[bp_indicies, ] %>%
    dplyr::mutate(side = rep(c("left", "right"), length(bp@cpts) - 1)) %>%
    dplyr::select(sample, seqnames, pos, side, smoothed_alt_freq)
}

#' Calculate Mode of a Numeric Vector
#'
#' Computes the mode of a given numeric vector, which is the value that
#'   appears most frequently.
#'
#' @param v Numeric vector.
#'
#' @return The mode of the vector.
get_mode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# alternative noise labelling
# add criteria that variant gets kept if it is in a run of length x or more
# (calculate RLE on the whole vector and then use index to figure out RLE group)

#' Determine if a Data Point is Noisy
#'
#' Checks whether a data point in a vector of allele frequencies is considered
#' noisy based on its neighborhood defined by a specified window size.
#' It evaluates both left and right windows separately and considers the point
#' noisy only if it does not match the mode of either window.
#'
#' @param freq_vector Numeric vector of allele frequencies.
#' @param index Integer, the index of the data point to check.
#' @param neighborhood_size Integer, the number of points on either side of
#'   the index to include in the neighborhood.
#'
#' @return Logical indicating whether the data point is noisy.
is_noisy <- function(freq_vector, index, neighborhood_size) {
  # Define the left and right window ranges, respecting the vector boundaries
  left_window_start <- max(1, index - neighborhood_size)
  left_window_end <- max(1, index - 1)
  right_window_start <- min(length(freq_vector), index + 1)
  right_window_end <- min(length(freq_vector), index + neighborhood_size)

  left_window = freq_vector[left_window_start:left_window_end]
  right_window = freq_vector[right_window_start:right_window_end]

  # Determine if the point matches either the left or right window's mode
  current_value <- freq_vector[index]

  # Calculate the mode of the left and right windows, if they are large enough,
  # and compare to current value
  left_window_match <- if (length(left_window) >= neighborhood_size)
    current_value == get_mode(left_window) else FALSE
  right_window_match <- if (length(right_window) >= neighborhood_size)
    current_value == get_mode(right_window) else FALSE

  is_noise <- !(left_window_match || right_window_match)

  return(is_noise)
}

#' Plot Breakpoints on Allele Frequency Data
#'
#' This function creates a plot of allele frequency data for each sample and chromosome,
#' with breakpoints overlaid as dashed vertical lines.
#'
#' @param allele_freq_df A data frame containing allele frequency data.
#' @param breakpoint_df A data frame containing breakpoint data.
#' @param allele_freq_stat Character string specifying the column name of the allele frequency
#'   statistic to plot. This should be an existing column in `allele_freq_df`.
#' @param sample_set Character vector specifying the samples to include in the plot.
#'   If empty, all samples will be included.
#' @param seqnames_set Character vector specifying the seqnames to include in the plot.
#'   If empty, all seqnames will be included.
#' @return A ggplot object displaying the allele frequency data with
#'   breakpoints overlaid, facilitating visual assessment of genomic variations.
#'
#' @importFrom ggplot2 ggplot aes geom_point coord_cartesian geom_vline facet_grid labs theme_minimal
#' @importFrom rlang sym
#'
#' @export
plot_breakpoints <- function(allele_freq_df,
                             breakpoint_df,
                             allele_freq_stat,
                             de_noise = FALSE,
                             sample_set = c(),
                             seqnames_set = c(),
                             point_size=0.7) {
  if (!allele_freq_stat %in% names(allele_freq_df)) {
    stop("allele_freq_stat must be a column name in allele_freq_df.")
  }

  if (length(sample_set) > 0) {
    allele_freq_df <- dplyr::filter(allele_freq_df, sample %in% sample_set)
    breakpoint_df <- dplyr::filter(breakpoint_df, sample %in% sample_set)
  }
  if (length(seqnames_set) > 0) {
    allele_freq_df <- dplyr::filter(allele_freq_df, seqnames %in% seqnames_set)
    breakpoint_df <- dplyr::filter(breakpoint_df, seqnames %in% seqnames_set)
  }
  if(de_noise){
    allele_freq_df <- dplyr::filter(allele_freq_df, !noisy_flag)
  }

  ggplot(allele_freq_df, aes(x = pos/1000, y = !!rlang::sym(allele_freq_stat),
                             color = cut(!!rlang::sym(allele_freq_stat),
                                         breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                                         labels = c("REF", "indeterminant", "HET", "indeterminant", "ALT")))) +
    geom_point(size=point_size) +
    scale_color_manual(values = c("REF" = "#1ABC9C",
                                  "HET" = "#F39C12",
                                  "ALT" = "#9B59B6",
                                  "indeterminant" = "black")) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_vline(data = breakpoint_df, aes(xintercept = pos/1000, group = side, color = side),
               linetype = 'dashed') +
    facet_grid(sample ~ seqnames, scales = 'free_x') +
    labs(y = "Allele Frequency", x = "Position (kb)", title = "Allele Frequency with Breakpoints", color = "Allele Frequency Stat")
}

plot_depth <- function(depth_gr,
                       stat = 'norm_depth',
                       sample_set = c(),
                       seqnames_set = c()) {

  depth_gr %>%
    as_tibble() %>%
    filter(sample %in% sample_set, seqnames %in% seqnames_set) %>%
    ggplot(aes(x = start /1000, xend = end /1000, y = norm_depth, yend = norm_depth, color = type)) +
    geom_segment(linewidth = 1) +  # Change the size parameter to adjust the thickness of the lines
    scale_color_manual(values = c('deletion' = 'blue', 'duplication' = 'red')) +
    coord_cartesian(ylim = c(0, 4)) +
    facet_grid(sample ~ seqnames, scales = 'free_x') +
    labs(y = stat, x = "Position (kb)", title = "CNVPytor Depth Stats", color = "CNV Type")
}


#all_samples_df = map(bsae@colData$sample, extract_alt_freq_dataframe) %>% bind_rows()

# for 1, 10, 15, 20
# [[1]]
# Time difference of 24.97729 secs
#
# [[2]]
# Time difference of 13.65013 secs
#
# [[3]]
# Time difference of 12.5191 secs
#
# [[4]]
# Time difference of 14.00582 secs

# x = map(c(10, 11, 12, 13, 14,15), ~{
#
#   # get start time
#   start_time = Sys.time()
#   map(bsae@colData$sample, extract_alt_freq_dataframe, ncpus=.)
#
#   # return end time
#   Sys.time() - start_time
#   })
#
# [[1]]
# Time difference of 11.17715 secs
#
# [[2]]
# Time difference of 14.55692 secs
#
# [[3]]
# Time difference of 13.6652 secs
#
# [[4]]
# Time difference of 13.26105 secs
#
# [[5]]
# Time difference of 13.80658 secs
#
# [[6]]
# Time difference of 13.26404 secs


extract_var_data = function(sample){
  # Get genotype data for sample C84_C84
  SeqArray::seqSetFilter(gds, sample.id = sample, variant.id = c8_df$variant_id)
  alleles <- SeqArray::seqGetData(gds, c(chr="chromosome",
                                         allele="allele",
                                         variant_id = "variant.id",
                                         depth='annotation/format/DP'))
  filter_df = tibble(seqnames = alleles$chr,
         variant_id = alleles$variant_id,
         allele = alleles$allele,
         depth = as.vector(alleles$depth)) %>%
    mutate(n_alt = str_count(allele, ',')) %>%
    filter(n_alt == 1,
           depth > 10,
           seqnames != 'CP022335.1') %>%
    separate(allele, c('ref', 'alt'), sep = ',') %>%
    filter(nchar(ref) == 1
           & nchar(alt) == 1)
  SeqArray::seqResetFilter(gds)

  SeqArray::seqSetFilter(gds, sample.id = sample,
                         variant.id = filter_df$variant_id)
  genotypes <- SeqArray::seqGetData(gds, "genotype")
  alleles <- SeqArray::seqGetData(gds, c(chr="chromosome",
                                         pos="position",
                                         allele="allele",
                                         variant_id = "variant.id",
                                         depth='annotation/format/DP',
                                         ad = 'annotation/format/AD'))
  ad_mat = matrix(alleles$ad$data, ncol = 2, byrow = TRUE)
  df = tibble(seqnames = alleles$chr,
         start = alleles$pos,
         variant_id = alleles$variant_id,
         allele = alleles$allele,
         depth = as.vector(alleles$depth)) %>%
    separate(allele, c('ref', 'alt'), sep = ',') %>%
    mutate(ad = ad_mat[,2],
           alt_freq = ad_mat[,2]/rowSums(ad_mat))
  SeqArray::seqResetFilter(gds)

  df %>%
    mutate(sample = sample)
}

classify_alt_freq = function(df, ref_freq_thres = 0.2,
                             hetero_freq_thres = c(0.4, 0.6),
                             alt_freq_thres = 0.8,
                             neighborhood_size = 20,
                             ...){
  # step 1: 'round' or smooth the allele frequency data to 0, 0.5 or 1 depending
  # on the thresholds
  message('step 1: classify variants to ref, het or alt')
  # TODO: remove filter, add flag. When calling is_noisy, do not consider
  #
  rounded_data <- df %>%
    mutate(pos = start) %>%
    dplyr::mutate(smoothed_alt_freq = case_when(
      alt_freq < ref_freq_thres ~ 0,
      alt_freq >= hetero_freq_thres[1] & alt_freq <= hetero_freq_thres[2] ~ 0.5,
      alt_freq > alt_freq_thres ~ 1,
      TRUE ~ NA_real_
    ))

  rounded_data_fltr = rounded_data %>%
    dplyr::filter(!is.na(smoothed_alt_freq))

  sample_chr_split = rounded_data_fltr %>%
    group_by(sample, seqnames) %>%
    group_split()

  # step 2: take the mode in the neighborhood around a given variant. If the
  # variant doesn't match the mode of a window of length neighborhood_size on
  # its left AND its right, it is labelled as noisy
  message('step 2: label noisy variants')
  doParallel::registerDoParallel(25)
  noisy_flag_list = foreach::foreach(
    df = sample_chr_split
  ) %dopar% {
    sapply(1:nrow(df), function(i) {
      is_noisy(df$smoothed_alt_freq, i, neighborhood_size)
    })
  }
  doParallel::stopImplicitCluster()

  rounded_data_fltr$noisy_flag = unlist(noisy_flag_list)

  rounded_data = rounded_data %>%
    left_join(select(rounded_data_fltr, seqnames, sample, pos, noisy_flag)) %>%
    replace_na(list(noisy_flag = TRUE))

  sample_chr_split2 = rounded_data %>%
    filter(!noisy_flag) %>%
    group_by(sample, seqnames) %>%
    group_split()

  # step 3: call breakpoints
  message('step 3: call breakpoints')
  # doParallel::registerDoParallel(1)
  # breakpoint_df_list = foreach::foreach(
  #   df = sample_chr_split2[54],
  #   .combine = list
  # ) %dopar% {
  #   de_noised_df = df %>% filter(!noisy_flag)
  #   if(nrow(de_noised_df) == 0){
  #     NULL
  #   } else{
  #     tryCatch({
  #       find_breakpoints(de_noised_df, ...)
  #     }, error = function(e){
  #       NULL
  #     })
  #   }
  #
  # }
  # doParallel::stopImplicitCluster()
  breakpoint_df_list = map(sample_chr_split2, find_breakpoints)

  list(allele_freq_df = rounded_data, breakpoint_df = bind_rows(breakpoint_df_list))

}

#' Temp Docs
#'
#' @param threshold_multiplier The number of median absolute deviations from the
#'  median to use as a threshold for filtering outlier regions. Default is 5.2
hampelFilter = function(g_prime_vector, hampel_threshold_multiplier = 5.2){
  ln_g_prime_vector <- log(g_prime_vector)

  median_ln_g_prime_vector <- median(ln_g_prime_vector)

  # calculate left median absolute deviation for the trimmed G' prime set
  MAD <-
    median(median_ln_g_prime_vector -
             ln_g_prime_vector[ln_g_prime_vector <= median_ln_g_prime_vector])

  # return a vector where the entry is TRUE if the original value IS AN OUTLIER
  (ln_g_prime_vector - median(ln_g_prime_vector)) > (hampel_threshold_multiplier * MAD)
}

# this is an alternate implementation -- need some research on
# whether to use this or to continue using current code
hampelFilter_alternate = function(g_prime_vector, hampel_threshold_multiplier = 5.2) {
  ln_g_prime_vector <- log(g_prime_vector)
  median_ln_g_prime_vector <- median(ln_g_prime_vector)

  # Calculate Median Absolute Deviation (MAD)
  MAD <- median(abs(ln_g_prime_vector - median_ln_g_prime_vector))

  threshold <- hampel_threshold_multiplier * MAD

  # return a vector where the entry is TRUE if the original value IS AN OUTLIER
  abs(ln_g_prime_vector - median_ln_g_prime_vector) > threshold
}

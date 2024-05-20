# library("depmixS4")
# library(here)
# library(dplyr)
#
# data("speed")
# set.seed(1)
#
# mod <- depmix(response = rt ~ 1,
#               data = speed,
#               nstates = 2,
#               trstart = runif(4))
#
# fm <- fit(mod)
#
# allele_freq_df = read_csv(here('temp/allele_freq_df_run_24.csv'))
#
# th0051_chr1 = allele_freq_df %>%
#   filter(sample == 'TH0051', seqnames == 'CP022321.1')
#
# mod <- depmix(alt_freq ~ 1,
#               data = th0051_chr1,
#               nstates = 2,
#               trstart = runif(4))
#
# fm <- fit(mod)
#
# th0051_chr1 %>%
#   mutate(postereriors = posterior(fm, type = 'local'))  %>%
#   mutate(posteriors = ifelse(posteriors == 2, 0, posteriors)) %>%
#   pivot_longer(cols = c('smoothed_alt_freq', 'posteriors'), names_to = 'metric', values_to = 'value')  %>%
#   ggplot(aes(pos, value)) +
#   geom_point() +
#   facet_grid(~metric)
#
# th0051_chr1 %>%
#   ggplot(aes(pos, smoothed_alt_freq)) +
#   geom_point()
#
#
# # Define a function to calculate the run length
# calculate_run_length <- function(values) {
#   rle_vals <- rle(values)$lengths
#   run_length <- rep(rle_vals, rle_vals)
#   return(run_length)
# }
#
# # Define a function to calculate the prior probability based on the run length
# calculate_prior <- function(run_length, decay_rate) {
#   p_change <- exp(-decay_rate * run_length)
#   return(p_change)
# }
#
# # Use dplyr to add the run length and prior probability to your tibble
# th0051_chr1 <- th0051_chr1 %>%
#   arrange(pos) %>% # Make sure the data is sorted by position
#   mutate(run_length = calculate_run_length(smoothed_alt_freq),
#          prior = calculate_prior(run_length, decay_rate = 0.001)) # Adjust decay_rate as needed
#
# mod <- depmix(smoothed_alt_freq ~ 1,
#               data = th0051_chr1,
#               nstates = 2,
#               family = gaussian(),
#               transition = ~ scale(prior),
#               instart = runif(4))
#
# fm <- fit(mod, emc=em.control(rand=FALSE))
#
# th0051_chr1 %>%
#   mutate(postereriors = posterior(fm, type = 'local'))  %>%
#   mutate(posteriors = ifelse(posteriors == 2, 0, posteriors)) %>%
#   pivot_longer(cols = c('smoothed_alt_freq', 'posteriors'), names_to = 'metric', values_to = 'value')  %>%
#   ggplot(aes(pos, value)) +
#   geom_point() +
#   facet_grid(~metric)



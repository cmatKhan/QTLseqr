# assays$alt_percentage <- assays$AD / assays$DP
# # create a boolean matrix where the value is 1 if the depth is greater than
# # 10 and the alt_percentage is greater than .9 or less than .1
# # if high_confidence_pl and/or high_confidence_gq are provided, use those
# # to label high confidence calls, also
# assays$high_confidence <- (assays$DP >= high_confidence_depth) &
#   ((assays$alt_percentage >= high_confidence_alt_percentage) |
#      (assays$alt_percentage <= (1 - high_confidence_alt_percentage)))
# if(is.numeric(high_confidence_pl)){
#   assays$high_confidence <- assays$high_confidence & (assays$PL <= high_confidence_pl)
# }
# if(is.numeric(high_confidence_gq)){
#   assays$high_confidence <- assays$high_confidence & (assays$GQ >= high_confidence_gq)
# }
# # replace any NA values in high_confidence with FALSE
# assays$high_confidence[is.na(assays$high_confidence)] <- FALSE

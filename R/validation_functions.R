
#' Validate that a given variable is a positive integer
#'
#' @param n The variable to validate
#' @param varname The name of the variable to include in the error message
#'
#' @keywords internal
validatePositiveInteger <- function(n, varname = "Variable") {
  if (!is.numeric(n) || n < 1 || floor(n) != n) {
    stop(sprintf("'%s' must be a positive integer", varname))
  }
}

# # validate .read_in_data arguments ----
# if(!is.logical(keep_multiallelic)){
#   error_messages <- c(error_messages,
#                       "`keep_multiallelic` must be a boolean")
# }


  # # validate assay related arguments ----
  # if(!is.numeric(high_confidence_depth) | high_confidence_depth < 0) {
  #   error_messages <- c(error_messages,
  #                       "`high_confidence_depth` must be a positive integer")
  # }
  # if(!is.numeric(high_confidence_alt_percentage) | high_confidence_alt_percentage < 0 | high_confidence_alt_percentage > 1) {
  #   error_messages <- c(error_messages,
  #                       "`high_confidence_alt_percentage` must be a number between 0 and 1")
  # }
  # if(!is.null(high_confidence_pl) && (!is.numeric(high_confidence_pl) | high_confidence_pl < 0)) {
  #   error_messages <- c(error_messages,
  #                       "`high_confidence_pl` must be a positive number or NULL")
  # }
  # if(!is.null(high_confidence_gq) && (!is.numeric(high_confidence_gq) | high_confidence_gq < 0)) {
  #   error_messages <- c(error_messages,
  #                       "`high_confidence_gq` must be a positive number or NULL")
  # }


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

#' @keywords internal
"_PACKAGE"

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Utility function to filter arguments
#' @param FUN a function
#' @param ... arguments to extract from the ... args based on FUN
#'
#' @return a list of arguments that are valid for FUN
filter_args <- function(FUN, ...) {
  args <- list(...)
  formals <- names(formals(FUN))
  args <- args[names(args) %in% formals]

  if(length(args) == 0){
    NULL
  } else{
    args
  }
}

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

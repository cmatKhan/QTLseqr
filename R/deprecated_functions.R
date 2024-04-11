#' Parse GT Table
#'
#' This function takes a data frame of GT data and a reference column,
#' compares each GT value to the reference, and returns a numeric matrix
#' where each element is 0 if the corresponding GT value is equal to the reference
#' and 1 otherwise.
#'
#' @param gt_data A data frame of GT data.
#' @param ref A character vector of reference values.
#' @param chunk_size An integer which sets the
#'
#' @return A numeric matrix where each element is 0 if the corresponding GT value
#' is equal to the reference and 1 otherwise.
#'
#' @importFrom purrr map
# @importFrom dplyr bind_cols
#'
#' @examples
#' gt_data <- data.frame(sample1 = c("A", "G", "C", "T", "A"),
#'                       sample2 = c("G", "A", "C", "T", "G"))
#' ref <- c("A", "G", "C", "T", "A")
#' .parse_gt_table(gt_data, ref)
#'
#' @keywords internal
.parse_gt_table <- function(gt_data, ref){
  # Iterate over columns in gt_data
  result_list <- purrr:::map(colnames(gt_data), function(col_name) {
    ifelse(gt_data[, col_name] == "./.", NA_integer_,
           ifelse(gt_data[, col_name] == ref, 0L, 1L))
  })

  # Store the result as a list, then use dplyr::bind_cols() %>% as.matrix
  names(result_list) = colnames(gt_data)
  dplyr::bind_cols(result_list) %>% as.matrix()
}

# test_that("parse_gt_table returns correct output", {
#   # Define the input data frame and reference vector
#   gt_data <- data.frame(sample1 = c("A", "G", "C", "T", "A"),
#                         sample2 = c("G", "A", "C", "T", "G"))
#   ref <- c("A", "G", "C", "T", "A")
#
#   # Call the function
#   result <- .parse_gt_table(gt_data, ref)
#
#   # Check that the result is a matrix with the correct dimensions
#   expect_true(is.matrix(result))
#   expect_equal(dim(result), c(length(ref), ncol(gt_data)))
#
#   # Check that the result has the correct values
#   expected_result <- matrix(c(0, 0, 0, 0, 0, 1, 1, 0, 0, 1), nrow = length(ref))
#   colnames(expected_result) = colnames(gt_data)
#   expect_equal(result, expected_result)
# })

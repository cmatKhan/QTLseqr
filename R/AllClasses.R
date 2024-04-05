#' @export
setClass(
  "BSAExperiment",
  contains = "RangedSummarizedExperiment",
  representation = representation(
    comparisons = "DataFrame"
  )
)

#' @import SummarizedExperiment
BSAExperiment <- function(comparisons = DataFrame(), ...) {

  # Create the SummarizedExperiment object
  se <- SummarizedExperiment(...)

  # Now create the BSAExperiment object, setting comparisons directly
  new("BSAExperiment", se, comparisons = comparisons)
}


S4Vectors::setValidity2("BSAExperiment", function(object) {
  msg <- NULL

  # check that assays has at least DP and AD
  if(is.null(object@assays)){
    warning("`assays` is null", call. = FALSE)
  } else if(!all(c("DP", "AD") %in% names(object@assays))) {
    msg <- c(msg, "The assays must contain at minimum 'DP' and 'AD'")
  }

  # Validate 'comparisons' DataFrame slot
  if (!inherits(object@comparisons, "DataFrame")) {
    msg <- c(msg, "'comparisons' slot must be a DataFrame.")
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

setClass(
  "BSAResults",
  contains = "DFrame",
  slots = c(
    rowRanges = "GRanges",
    analysisParams = "list"
  )
)

BSAResults <- function(DataFrame = S4Vectors::DataFrame(),
                          rowRanges = GenomicRanges::GRanges(),
                          analysisParams = list()) {
  # Create and return the BSAResults object with rowData
  new("BSAResults",
      DataFrame,
      rowRanges = rowRanges,
      analysisParams = analysisParams
  )
}

S4Vectors::setValidity2("BSAResults", function(object) {
  msg <- NULL

  # Validate 'rowData' slot
  if (!inherits(object@rowRanges, "GRanges")) {
    msg <- c(msg, "'rowData' slot must be a GRanges object.")
  }

  # Validate 'analysisParams' slot
  if (!is.list(object@analysisParams)) {
    msg <- c(msg, "'analysisParams' slot must be a list.")
  }

  # Validate 'DFrame' slot
  if (!inherits(object, "DataFrame")) {
    msg <- c(msg, "'DataFrame' slot must be a DataFrame.")
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

# Define the constructor function for the S3 object
createPopulationDepthList <- function(depth_list) {
  # Validate that the input is a list
  if (!is.list(depth_list)) {
    stop("Input is not a list")
  }

  required_top_names = c('population_1', 'population_2')
  if(!all(required_top_names %in% names(depth_list))) {
    stop('The first level of names must be `population_1` and `population_2`')
  }

  # Validate that the list has the correct names
  required_population_names <- c("ref",  "alt")
  for(population in depth_list){
    if(!all(required_population_names %in% names(population))) {
      stop(population, " list must have names `ref` and `alt`.")
    }
  }

  # Validate that all of the elements are the same length
  n = length(depth_list$population_1$ref)
  if(!all(sapply(unlist(depth_list, recursive=FALSE), length) == n)) {
    stop("Input depth vectors for both populations ",
         "ref/alt must all be the same length. At least 1 is not.")
  }

 # Validate that all of the elements are numeric
  if (!all(sapply(unlist(depth_list), is.numeric))) {
    stop("Input list elements are not numeric")
  }

  # If validation passes, set the class attribute and return the object
  class(depth_list) <- "PopulationDepthList"

  depth_list
}


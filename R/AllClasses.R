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

  # Capture all additional arguments in a list
  argsList <- list(...)

  # Check if 'metadata' is explicitly passed
  if("metadata" %in% names(argsList)) {
    # Access the metadata
    metadata <- argsList$metadata
  } else {
    metadata <- list()
  }

  if(!is.list(metadata)){
    warning('metadata was either not passed explicitly or is not a list. ',
            'Setting to an empty list with null initial values for ',
            '`population_1_n`, `population_2_n` and `population_structure`')
    metadata <- list()
  }

  if(is.null(metadata$population_1_n)){
    metadata$population_1_n <- NA
  }

  # if population_1_n exists, but population_2_n is null, then set population_2_n
  # to population_1_n value
  if (!is.null(metadata$population_1_n) && is.null(metadata$population_2_n)) {
    warning("setting `population_2_n` to `population_1_n` value")
    metadata$population_2_n <- metadata$population_1_n
  }

  # if `population_structure` doesn't exist, then initialize it to NULL
  if(is.null(metadata$population_structure)){
    metadata$population_structure <- NA
  }

  # Put the modified metadata back into argsList
  argsList$metadata <- metadata

  # Create the SummarizedExperiment object
  se <- do.call("SummarizedExperiment", argsList)

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

  if(rlang::is_empty(object@metadata)){
    msg <- c(msg, paste0("The metadata slot is empty. `population_1_n`, ",
             "`population_2_n` and `population_structure` ",
             "must be set in the metadata"))
  }

  # validate that the names 'population_1_n', 'population_2_n' and 'population_structure'
  # are in object@metadata
  if(!all(c("population_1_n", "population_2_n", "population_structure") %in% names(object@metadata))){
    msg <- c(msg, paste0("The metadata is missing the following fields: ",
             setdiff(c("population_1_n", "population_2_n", "population_structure"),
                     names(object@metadata))))
  }

  # if population_1_n and population_2_n are not NA, then check their values
  for (field in c("population_1_n", "population_2_n")) {
    if (!is.na(object@metadata[[field]])) {
      tryCatch({
        validatePositiveInteger(object@metadata[[field]], field)
      }, error = function(e) {
        if (is.na(object@metadata[[field]])) {
          warning(paste0("`", field, "` is `NA`. This must be set to a valid positive integer ",
                         "prior to running certain operations"))
        } else {
          msg <- c(msg, e$message)
        }
      })
    }
  }

  # if population_structure is not NA, check the value
  if(!is.na(object@metadata$population_structure)){
    if(!object@metadata$population_structure %in% c("F2", "RIL")){
      msg <- c(msg, paste0("The metadata slot must contain ",
                           "'population_structure' as either 'F2' or 'RIL'. ",
                           "Other population structures not currently supported"))
    }
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

#' Population Depth List Object Constructor
#' @export
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


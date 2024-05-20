# hard coded parameters, eg colnames, used in multiple methods ----
methods_params <- list(
  alt_depth_slot = "AD",
  total_depth_slot = "DP",
  variant_filter_col = "variant_filter",
  sample_filter_col = "sample_filter",
  matrix_filter_slot = "entry_filter"
)

# Local Generics ----

## Common params ----

#' @name commonParams
#' @keywords internal
#' @param object A `BSAExperiment` object
#' @param population_1_sample The name of the sample used as the first
#'   population. In the context of the Takagi 2013 paper, this would be the
#'   'low bulk' sample.
#' @param population_2_sample The name of the sample used for the second
#'   population. In the context of the Takagi 2013 paper, this would be the
#'   'high bulk' sample.
#' @param sample_set A character vector of samples, possibly of length 1.
#'   All entries in `sample_set` must be present in the `sample` column of
#'   the colData
NULL

## comparisons ----

#' Getter/Setter of the BSAExperiment comparison slot
#'
#' @inheritParams commonParams
#' @param value A dataframe describing comparisons between two populations.
#'   Must have two columns where the entries in each row are names of samples
#'   in the `sample` column of the colData
#'
#' @return The setter will return a BSAExperiment object with an updated
#'   comparison slot while the getter will return the comparison slot's
#'   DataFrame
#'
#' @rdname comparisons
#'
#' @export
setReplaceMethod("comparisons", "BSAExperiment", function(object, value) {
  # Update object
  object@comparisons <- value
  # validate
  validObject(object)
  # return
  object
})

#' @rdname comparisons
#' @export
setMethod("comparisons", "BSAExperiment", function(object) {
  object@comparisons
})

## totalDepth ----

#' Getter/Setter of the BSAExperiment totalDepth assay matrix
#'
#' Set the totalDepth matrix in the BSAExperiment object, or retrieve it.
#'   For the getter, a vector of BSAEXperiment colnames, possibly of length 1,
#'   may be passed to extract a subset of the samples from the BSAExperiment
#'   object
#'
#' @inheritParams commonParams
#' @param value A numeric matrix storing the total depth at each variant
#'   location in the BSAExperiment. The dimensions of this matrix must be
#'   equal to the dimension of the BSAExperiment, and must have the same
#'   row and column names
#'
#' @return The setter will return a BSAExperiment object with the
#'   total depth slot updated. The getter will return a matrix of
#'   total depths
#'
#' @rdname totalDepth
#'
#' @export
setReplaceMethod("totalDepth", "BSAExperiment", function(object, value) {
  # Update object
  tryCatch(
    {
      assay(object, methods_params$total_depth_slot) <- value
    },
    error = function(e) {
      stop("Error adding Total Depth data: ", e$message)
    }
  )
  # validate
  validObject(object)
  # return
  object
})

#' @inheritParams commonParams
#' @rdname totalDepth
#' @export
setMethod("totalDepth", "BSAExperiment", function(object, sample_set = NULL) {
  if (is.null(sample_set)) {
    local_sample_set <- colnames(object)
  } else {
    if (!all(sample_set %in% colnames(object))) {
      stop(
        "At least one of the provided samples is not present in the colnames ",
        "of the BSAExperiment object"
      )
    }
    local_sample_set <- sample_set
  }
  assay(object[, local_sample_set], methods_params$total_depth_slot)
})

## altDepth ----

#' Getter/Setter of the BSAExperiment altDepth assay matrix
#'
#' Set the altDepth matrix in the BSAExperiment object, or retrieve it.
#'   For the getter, a vector of BSAEXperiment colnames, possibly of length 1,
#'   may be passed to extract a subset of the samples from the BSAExperiment
#'   object
#'
#' @inheritParams commonParams
#' @param value A numeric matrix storing the alt depth at each variant
#'   location in the BSAExperiment. The dimensions of this matrix must be
#'   equal to the dimension of the BSAExperiment, and must have the same
#'   row and column names
#'
#' @return The setter will return a BSAExperiment object with the
#'   alt depth slot updated. The getter will return a matrix of
#'   alt depths
#'
#' @rdname altDepth
#'
#' @export
setReplaceMethod("altDepth", "BSAExperiment", function(object, value) {
  # Update object
  tryCatch(
    {
      assay(object, methods_params$alt_depth_slot) <- value
    },
    error = function(e) {
      stop("Error adding Alt Depth data: ", e$message)
    }
  )
  # validate
  validObject(object)
  # return
  object
})

#' @inheritParams commonParams
#' @rdname altDepth
#' @export
setMethod("altDepth", "BSAExperiment", function(object, sample_set = NULL) {
  if (is.null(sample_set)) {
    local_sample_set <- colnames(object)
  } else {
    if (!all(sample_set %in% colnames(object))) {
      stop(
        "At least one of the provided samples is not present in the colnames ",
        "of the BSAExperiment object"
      )
    }
    local_sample_set <- sample_set
  }
  assay(object[, local_sample_set], methods_params$alt_depth_slot)
})

## refDepth ----

#' Get the Reference Depths from the BSAExperiment
#'
#' Optionally, a vector of BSAEXperiment colnames, possibly of length 1,
#'   may be passed to extract a subset of the samples from the BSAExperiment
#'   object
#'
#' @inheritParams commonParams
#'
#' @return A matrix of total depths
#'
#' @rdname refDepth
#' @export
setMethod("refDepth", "BSAExperiment", function(object, sample_set = NULL) {
  tryCatch({
    totalDepth(object, sample_set) - altDepth(object, sample_set)
    }, error = function(e){
      stop(e$message)
    })
})

## variantFilter ----

#' Getter/Setter of the BSAExperiment variantFilter vector
#'
#' The variantFilter is a boolean vector which is added as a column to the
#'   rowRanges mcols. A value of `FALSE` indicates that that variant should
#'   be removed from the filtered dataset. When passing a vector to the setter,
#'   it must be equal to the number of rows in the BSAExperiment object
#'
#'
#' @inheritParams commonParams
#' @param value A boolean vector of length equal to the rows of the
#'   BSAExperiment object. A value of `FALSE` indicates that that variant
#'   should be removed from the filtered dataset.
#'
#' @return The setter will return a BSAExperiment object with the variant
#'   filter vector added to the mcols of the rowRanges. The getter will return
#'   the boolean vector
#'
#' @rdname variantFilter
#'
#' @export
setReplaceMethod("variantFilter", "BSAExperiment", function(object, value) {
  if (!is.logical(value)) {
    stop("The value must be a logical vector.")
  }
  if (length(value) != nrow(object)) {
    stop(
      "The length of the value must be equal to the number of ",
      "rows in the BSAExperiment object."
    )
  }

  # Update object
  rowData(object)[[methods_params$variant_filter_col]] <- value
  # Validate object
  validObject(object)
  # Return updated object
  object
})


#' @inheritParams commonParams
#' @rdname variantFilter
#' @export
setMethod("variantFilter", "BSAExperiment", function(object) {
  if (!methods_params$variant_filter_col %in% colnames(rowData(object))) {
    message(
      "There is no explicit filter column in the rowData. ",
      "Returning TRUE for all variants."
    )
    return(rep(TRUE, nrow(object)))
  } else {
    return(rowData(object)[[methods_params$variant_filter_col]])
  }
})

## sampleFilter ----

#' Getter and Setter of the BSAExperiment sampleFilter vector
#'
#' The sampleFilter is a boolean vector which is added as a column to the
#'   colData. A value of `FALSE` indicates that that entire sample should
#'   be removed from the filtered dataset. When passing a vector to the setter,
#'   it must be equal to the number of columns in the BSAExperiment object
#'
#' @inheritParams commonParams
#' @param value A boolean vector of length equal to the columns of the
#'   BSAExperiment object. A value of `FALSE` indicates that that sample
#'   should be removed from the filtered dataset.
#'
#' @return The setter will return a BSAExperiment object with the sample
#'   filter vector added to the colData. The getter will return
#'   the boolean vector
#'
#' @rdname sampleFilter
#'
#' @export
setReplaceMethod("sampleFilter", "BSAExperiment", function(object, value) {
  if (!is.logical(value)) {
    stop("The value must be a logical vector.")
  }
  if (length(value) != ncol(object)) {
    stop(
      "The length of the value must be equal to the number of ",
      "columns in the BSAExperiment object."
    )
  }

  # Update object
  colData(object)[[methods_params$sample_filter_col]] <- value
  # Validate object
  validObject(object)
  # Return updated object
  object
})

#' @inheritParams commonParams
#' @rdname sampleFilter
#' @export
setMethod("sampleFilter", "BSAExperiment", function(object) {
  if (!methods_params$sample_filter_col %in% colnames(colData(object))) {
    message(
      "There is no explicit filter column in the colData. ",
      "Returning TRUE for all samples."
    )
    return(rep(TRUE, ncol(object)))
  } else {
    return(colData(object)[[methods_params$sample_filter_col]])
  }
})

## matrixFilter ----

#' Getter and Setter of the BSAExperiment matrixFilter matrix
#'
#' The matrixFilter is a boolean matrix of `nrow(object) x ncol(object)` that
#'   is stored in the assays slot. A value of `FALSE` indicates that the
#'   the entry in the other matricies in the assays slot that correspond to
#'   the entry in the matrixFilter should be set to `NA`.
#'
#' @inheritParams commonParams
#' @param value A boolean matrix of dim .
#'   A value of `FALSE` indicates that the corresponding entry of all of the
#'   assay data matrices should be set to NA.
#'
#' @return The setter will return a BSAExperiment object with the matrix
#'   filter matrix added to the assays The getter will return
#'   the boolean matrix
#'
#' @rdname matrixFilter
#' @export
setReplaceMethod("matrixFilter", "BSAExperiment", function(object, value) {
  if (!is.matrix(value)) {
    stop("The value must be a matrix.")
  }
  if (!is.logical(value)) {
    stop("The value must be a logical matrix.")
  }
  if (!all(dim(value) == c(nrow(object), ncol(object)))) {
    stop(
      "The dimensions of the value must match the dimensions of ",
      "the BSAExperiment object."
    )
  }

  # Update object
  assays(object)[[methods_params$matrix_filter_slot]] <- value
  # Validate object
  validObject(object)
  # Return updated object
  object
})

#' @inheritParams commonParams
#' @rdname matrixFilter
#' @export
setMethod("matrixFilter", "BSAExperiment", function(object) {
  if (!methods_params$matrix_filter_slot %in% names(assays(object))) {
    message(
      "There is no explicit filter matrix in `assays`. ",
      "Returning a matrix of FALSE."
    )
    matrix(FALSE, nrow(object), ncol(object))
  } else {
    assay(object, methods_params$matrix_filter_slot)
  }
})

## applyFilter ----

#' Apply the variant, sample and/or matrix filters to a BSAExperiment object
#'
#' If any one of the selected, or by default, all, of the filters are set,
#'   the corresponding entries in the BSAExperiment object will either be
#'   removed, in the case of the variant or sample filter, or all entries
#'   in the corresponding assays slots will be set to NA.
#'
#' @inheritParams commonParams
#' @param variant_filter Boolean, default `TRUE`. Set to FALSE to use the
#'   variant filter to remove variants from the BSAExperiment object.
#' @param sample_filter Boolean, default `TRUE`. Set to FALSE to use the
#'   sample filter to remove samples from the BSAExperiment object.
#' @param matrix_filter Boolean, default `TRUE`. Set to FALSE to use the
#'   matrix filter to remove entries from the BSAExperiment object.
#'
#' @return A BSAExperiment object filtered based on the supplied variant,
#'   sample and/or entry filters
#'
#' @rdname applyFilter
#' @export
setMethod("applyFilter", "BSAExperiment", function(object,
                                                    variant_filter = TRUE,
                                                    sample_filter = TRUE,
                                                    matrix_filter = TRUE) {
  local_object = object

  # matrix filter must be applied first because the variant and sample filters
  # may change the dimension of the object. Where the entry filter is FALSE,
  # values in the assay data will be replaced with NA
  if (matrix_filter) {
    entry_filter <- matrixFilter(local_object)
    data_slots <- setdiff(names(assays(local_object)),
                          methods_params$entry_filter_slot)
    for (slot in data_slots) {
      assay_data <- assay(local_object, slot)
      assay_data[!entry_filter] <- NA
      assay(local_object, slot) <- assay_data
    }
  }

  if (variant_filter){
    local_object = local_object[variantFilter(local_object), ]
  }

  if (sample_filter){
    local_object = local_object[, sampleFilter(local_object)]
  }

  local_object
})

## populationDepths ----

#' Extract population alt and ref allele depths from a BSAExperiment object
#'
#' @inheritParams commonParams
#'
#' @rdname populationDepths
#'
#' @export
setMethod(
  "populationDepths", "BSAExperiment",
  function(object, population_1_sample, population_2_sample) {
    # validate that population_1_sample is in the colData(object)$sample
    for (sample in c(population_1_sample, population_2_sample)) {
      if (!sample %in% colnames(object)) {
        stop(paste0("sample: ", sample, " not found in colnames(object)"))
      }
    }

    # create a list of population ALT depths
    tmp <- list(
      population_1 = list(
        alt = assays(object[, population_1_sample])$AD
      ),
      population_2 = list(
        alt = assays(object[, population_2_sample])$AD
      )
    )

    # add the REF depths to each of the populations
    tmp$population_1$ref <- assays(object[, population_1_sample])$DP - tmp$population_1$alt

    tmp$population_2$ref <- assays(object[, population_2_sample])$DP - tmp$population_2$alt

    createPopulationDepthList(tmp)
  }
)

## snpIndex -----

#' Calculate the SNP index of a given sample
#'
#' @inheritParams commonParams
#'
#' @rdname snpIndex
#'
#' @export
setMethod(
  "snpIndex", "BSAExperiment",
  function(object, sample_set = NULL) {
    tryCatch({
      altDepth(object, sample_set) / totalDepth(object, sample_set)
    }, error = function(e){
      stop(e$message)
    })
  }
)

## deltaSnpIndex ----

#' Calculate \eqn{\Delta}{Delta}(Alt Allele Frequency) Between Two Populations
#'
#' Calculate the difference in alternate allele frequency between population_2
#' and population_1
#'
#' @inheritParams commonParams
#'
#' @rdname deltaSnpIndex
#'
#' @return A numeric vector that results from
#' \eqn{Population_2\ Alternate\ Allele\ Frequency\ -\ Population_1\ Alternate\ Allele\ Frequency}
#'
#' @examples
#' # Assuming `bsae` is a BSAExperiment object with populations "Pop1"
#' # and "Pop2"
#' delta_freqs <- deltaSnpIndex(bsae, "Pop1", "Pop2")
#'
#' # View the first few differences in allele frequency
#' head(delta_freqs)
#'
#' @export
setMethod(
  "deltaSnpIndex", "BSAExperiment",
  function(object, population_1_sample, population_2_sample) {
    if(length(population_1_sample) != length(population_2_sample)){
      stop("The number of samples in population_1 ",
           "and population_2 must be equal")
    }
    # Retrieve depth information for both populations
    population_1_snp_index <- snpIndex(object, population_1_sample)

    population_2_snp_index <- snpIndex(object, population_2_sample)

    # Compute the difference in allele frequencies between population_2
    # and population_1
    delta_snp_index <- population_2_snp_index - population_1_snp_index

    # Label the columns of the delta frequency dataframe to reflect
    # the comparison
    colnames(delta_snp_index) <- paste(population_2_sample,
      population_1_sample,
      sep = "-vs-"
    )

    # Return the dataframe containing delta allele frequencies
    delta_snp_index
  }
)

# Outside Generics ----

## summary ----

#' Summarize a BSAResults object
#'
#' @inheritParams commonParams
#' @rdname summary
#'
#' @export
setMethod("summary", "BSAResults", function(object) {
  message("Summary of BSAResults Object")

  # Summary of rowData (GRanges)
  if (is(object@rowData, "GRanges")) {
    message("\nrowData (GRanges):")
    gr <- object@rowData
    cat("Columns: ", names(mcols(gr)), "\n")
    cat("Ranges info: seqnames, ranges, strand (if present)\n")
  }

  # Numeric summaries of DataFrame columns
  message("\nDataFrame Summary:")
  df <- as.data.frame(object)
  summaryCols <- sapply(df, is.numeric)
  summaryDf <- summary(df[, summaryCols, drop = FALSE])
  print(summaryDf)

  # Analysis parameters
  message("\nAnalysis Parameters:")
  print(object@analysisParams)
})

## aggregate ----

#' Aggregate Data in BSAExperiment
#'
#' Aggregates assay data in a BSAExperiment object based on specified
#'   grouping factors in colData. The aggregation is performed by a
#'   user-specified function, typically a sum, mean, etc.
#'
#' @param x An object of class BSAExperiment.
#' @param by A formula indicating the grouping factors,
#'   e.g., ~ condition + sac_day.
#' @param FUN Function to apply to each group for aggregation.
#'
#' @importFrom stats aggregate
#' @importFrom rlang is_formula
#' @importFrom S4Vectors SimpleList
#'
#' @rdname aggregate
#'
#' @export
setMethod("aggregate", "BSAExperiment", function(x, by, FUN, ...) {
  # Check that `by` is a formula
  if (!rlang::is_formula(by)) {
    stop("`by` must be a formula.")
  }

  # Extract the variables from the formula and check their presence and
  # type in colData
  vars <- all.vars(by)
  if (!all(vars %in% colnames(colData(x)))) {
    missing_vars <- vars[!vars %in% colnames(colData(x))]
    stop(
      "The following variables in 'by' are not in colData: ",
      paste(missing_vars, collapse = ", "), "."
    )
  }

  # check that FUN is a function
  if (!is.function(FUN)) {
    stop("`FUN` must be a function, eg rowSums or rowMeans.")
  }

  # Step 1: Create a model frame to store group information according to
  # the `by` formula
  model_frame <- model.frame(by, data = colData(x))
  interaction_vector <- droplevels(interaction(model_frame))

  # Step 2: Aggregate assay data
  assay_matrix_names <- c("AD", "DP")
  agg_assays <- lapply(assay_matrix_names, function(n) {
    assay_matrix <- x@assays@data[[n]]

    agg_matrix <- matrix(
      nrow = nrow(assay_matrix),
      ncol = length(levels(interaction_vector))
    )

    rownames(agg_matrix) <- rownames(assay_matrix)
    colnames(agg_matrix) <- levels(interaction_vector)

    filtered_args <- filter_args(FUN, ...)

    for (i in levels(interaction_vector)) {
      samples_to_combine <- colnames(assay_matrix)[interaction_vector == i]
      # Apply FUN with filtered arguments
      agg_matrix[, i] <- do.call(FUN, c(list(assay_matrix[, samples_to_combine, drop = FALSE]), filtered_args))
    }
    agg_matrix
  })

  # name the vectors and cast to S4Vectors::SimpleList
  names(agg_assays) <- assay_matrix_names
  agg_assays <- S4Vectors::SimpleList(agg_assays)

  agg_coldata <- DataFrame(
    condition = levels(interaction_vector)
  )
  rownames(agg_coldata) <- agg_coldata$condition

  new_bsae <- BSAExperiment(
    assays = agg_assays,
    rowRanges = rowRanges(x),
    colData = agg_coldata,
    metadata = x@metadata
  )

  message(
    "You will may want to re-do the `comparisons` frame and add ",
    "variables to the colData"
  )

  new_bsae
})

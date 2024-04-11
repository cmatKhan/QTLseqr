#' @name commonParams
#' @param object A `BSAExperiment` object
#' @param population_1_sample The name of the sample used as the first population.
#'   In the context of the Takagi 2013 paper, this would be the 'low bulk' sample.
#' @param population_2_sample The name of the sample used for the second population.
#'   In the context of the Takagi 2013 paper, this would be the 'high bulk' sample.
NULL


#' Get the BSAExperiment comparisons frame
#'
#' @inheritParams commonParams
#'
#' @rdname BSAExperiment-methods
#' @export
setMethod("comparisons", "BSAExperiment", function(object) {
  object@comparisons
})

#' Set the BSAExperiment comparison frame
#'
#' @inheritParams commonParams
#' @param value A dataframe describing comparisons between two populations.
#'   Must have two columns where the entries in each row are names of samples
#'   in the `sample` column of the colData
#'
#' @return a BSAExperiment object with the comparisons slot updated
#'
#' @rdname BSAExperiment-methods
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

#' Create Sample Comparison Frame
#'
#' Create a dataframe which describes comparisons between a given
#'   condition and all other conditions described in a pool dataframe
#'
#' @inheritParams commonParams
#'
#' @param grouping_variable column by which to group rows of the dataframe.
#'   Defaults to 'batch'.
#' @param population_variable column which describes the various conditions
#'   of each sample, eg if tissue, then the levels of this column might be
#'   c(lung, brain, YPD, inoculum), Default: 'condition'
#' @param population_base_condition condition against which to compare all
#'   other conditions, Default: 'inoculum'
#' @param var1_name the output frame will have two columns, the first storing
#'   the samples which correspond to the `population_base_condition`,
#'   the other to the other sample conditions, Default: 'population_1'
#' @param var2_name as in var1_name, this will rename the second column
#'   in the output frame, Default: 'population_2'
#' @param base_cond_in_each_group whether to include the base condition in
#' each group. Default TRUE
#'
#' @return A two column dataframe where the first column is the condition
#'   against which all other sample conditions in that group are compared.
#'   An example structure is:
#'   \tabular{rcc}{ \tab population_2 \tab population_1 \cr \tab \eqn{P1.1I}
#'   \tab \eqn{P1.1L} \cr \tab
#' \eqn{P1.1I} \tab \eqn{P1.1B} \cr}
#' @details This prepares samples for QTLseqr
#' @examples
#' if (interactive()) {
#'   library(dplyr)
#'   # NOTE! "P2.1I" is a singleton
#'   sample_example <- c(
#'     "P1.1I", "P1.1L", "P1.1Y",
#'     "P1.2B", "P1.2I", "P1.2L", "P1.2Y", "P2.1I"
#'   )
#'   pool_construction <- tibble(sample = sample_example) %>%
#'     mutate(group = str_remove(sample, "\\w$")) %>%
#'     mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}I"), "inoculum", NA)) %>%
#'     mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}Y"), "ypd", cond)) %>%
#'     mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}L"), "lung", cond)) %>%
#'     mutate(cond = ifelse(str_detect(sample, "P[[:alnum:]].{1,3}B"), "brain", cond)) %>%
#'     mutate(bulk = ifelse(cond == "inoculum", "low", "high"))
#'   sample_comparison_frame(pool_construction)
#' }
#'
#' @rdname BSAExperiment-methods
#'
#' @importFrom dplyr setdiff group_by group_split pull filter distinct mutate rename as_tibble
#' @importFrom rlang sym
#' @importFrom purrr map
#' @importFrom S4Vectors DataFrame
#'
#' @export
setMethod(
  "createComparisonsFrame",
  "BSAExperiment",
  function(object,
           grouping_variable = "group",
           population_variable = "bulk",
           population_base_condition = "low",
           var1_name = "population_1",
           var2_name = "population_2",
           base_cond_in_each_group = TRUE) {
    df <- colData(object) %>%
      dplyr::as_tibble()

    # check colnames of dataframe against expected colnames
    expected_colnames <- c(grouping_variable, population_variable, "sample")

    if (any(!expected_colnames %in% colnames(df)) & base_cond_in_each_group) {
      stop(paste0(
        "input dataframe colnames must possess a field named sample, ",
        "and at least the `grouping_variable` and `population_variable` ",
        sprintf(
          "in any order: %s",
          paste(expected_colnames, collapse = ",")
        )
      ))
    }

    grouping_variable_symbol <- rlang::sym(grouping_variable)
    population_variable_symbol <- rlang::sym(population_variable)

    # group the dataframe on grouping_variable, and then split into a list of
    # dataframes where each item is one group
    if (base_cond_in_each_group) {
      group_split_df <- df %>%
        dplyr::group_by(!!grouping_variable_symbol) %>%
        dplyr::group_split()

      df_list <- purrr::map(group_split_df, function(df) {
        comparison_condition_sample <- df %>%
          dplyr::filter(!!population_variable_symbol == population_base_condition) %>%
          dplyr::pull(sample)

        expand_grid_df <- expand.grid(comparison_condition_sample, df$sample) %>%
          dplyr::distinct(Var2, .keep_all = TRUE) %>%
          dplyr::mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
          dplyr::filter(Var1 != Var2)

        expand_grid_df
      })

      # combine the list of tables into a single table again
      result_df <- do.call("rbind", df_list) %>%
        # rename the Var1 and Var2
        dplyr::rename(
          !!rlang::sym(var1_name) := Var1,
          !!rlang::sym(var2_name) := Var2
        )
    } else {
      population_1_condition <- df %>%
        dplyr::filter(!!population_variable_symbol == population_base_condition) %>%
        dplyr::pull(sample)

      if (length(population_1_condition) > 1) {
        stop(
          "There is more than 1 low population condition -- ",
          "setting base_cond_in_each_group failed. Consider splitting up ",
          "the dataframe and running this function on parts, or set ",
          "base_cond_in_each_group to TRUE and allow those conditions which ",
          "do not have base conditions to be dropped."
        )
      }

      population_2_conditions <- df %>%
        dplyr::filter(!!population_variable_symbol !=
          population_base_condition) %>%
        dplyr::pull(sample)

      result_df <- expand.grid(population_1_condition, population_2_conditions) %>%
        dplyr::rename(
          !!rlang::sym(var1_name) := Var1,
          !!rlang::sym(var2_name) := Var2
        )
    }

    object@comparisons <- S4Vectors::DataFrame(result_df)

    object
  }
)

#' Summarize a BSAResults object
#'
#' @inheritParams commonParams
#' @rdname BSAExperiment-methods
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

#' Extract population alt and ref allele depths from a BSAExperiment object
#'
#' @inheritParams commonParams
#'
#' @rdname BSAExperiment-methods
#'
#' @export
setMethod(
  "populationDepths", "BSAExperiment",
  function(object, population_1_sample, population_2_sample) {
    # validate that population_1_sample is in the colData(object)$sample
    if (!(population_1_sample %in% colData(object)$sample)) {
      stop(paste0("population_1_sample: ", population_1_sample, " not found in colData(object)$sample"))
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

#' Calculate \eqn{\Delta}{Delta}(Alt Allele Frequency) Between Two Populations
#'
#' Calculate the difference in alternate allele frequency between population_2
#' and population_1
#'
#' @inheritParams commonParams
#'
#' @return A numeric vector that results from
#' \eqn{Population_2\ Alternate\ Allele\ Frequency\ -\ Population_1\ Alternate\ Allele\ Frequency}
#'
#' @examples
#' # Assuming `bsae` is a BSAExperiment object with populations "Pop1" and "Pop2"
#' delta_freqs <- deltaAltFrequency(bsae, "Pop1", "Pop2")
#'
#' # View the first few differences in allele frequency
#' head(delta_freqs)
#'
#' @export
setMethod(
  "deltaAltFrequency", "BSAExperiment",
  function(object, population_1_sample, population_2_sample) {
    # Retrieve depth information for both populations
    depths <- populationDepths(object, population_1_sample, population_2_sample)

    # Calculate allele frequency for population 1 by dividing ALT depth by total depth (DP)
    population_1_alt_frequency <- (depths$population_1$alt / assays(object[, population_1_sample])$DP)

    # Calculate allele frequency for population 2 similarly
    population_2_alt_frequency <- (depths$population_2$alt / assays(object[, population_2_sample])$DP)

    # Compute the difference in allele frequencies between population_2 and population_1
    delta_alt_frequency <- population_2_alt_frequency - population_1_alt_frequency

    # Label the columns of the delta frequency dataframe to reflect the comparison
    colnames(delta_alt_frequency) <- paste(colnames(depths$population_2$alt),
      colnames(depths$population_1$alt),
      sep = "-vs-"
    )

    # Return the dataframe containing delta allele frequencies
    delta_alt_frequency
  }
)

#' Create a BSAExperiment object from a GATK `varianttotable` Table
#'
#' This function reads a GATK `varianttotable` output table and a
#'   column data file to create a BSAExperiment object. The GATK
#'   `VariantsToTable` command should be used to generate the input table. See
#'   the GATK `varianttotable` documentation for details.
#'
#' @inheritParams .read_in_gatk_table
#'
#' @param col_data_path Optional. If null, samplenames from the gatk table will
#'   be used to create a single column of metadata with the column name `sample`.
#'   However, additional sample metdata can be provided with this argument. If
#'   not `NULL`, this must be a character string specifying the path to the
#'   sample metadata. This must be a csv file with the extention `.csv`. The
#'   column `sample` must exist, and it must correspond to the sample names
#'   in the GATK `varianttotable` table.
#' @param metadata A list containing metadata for the BSAExperiment object.
#'   Recognized entries are `population_structure`, which may be 'F2' for an
#'   F2 population and 'RIL' for Recombinant Inbred Lines, `population_1_n`, which
#'   are the number of individuals in the 'low bulk' sample, and `population_2_n`,
#'   which are the number of individuals in the 'high bulk' sample. If `population_2_n`
#'   is not provided, it is equal to `population_1_n` by default.
#'
#' @details
#' To create the input table for this function, use the following GATK `VariantsToTable` command:
#' \preformatted{
#'
#' gatk VariantsToTable \
#'    -V input.vcf \
#'    --split-multi-allelic \
#'    -F CHROM -F POS -F REF -F ALT -F MULTI-ALLELIC -F TYPE -F QUAL \
#'    -GF AD -GF DP -GF PL -GF GQ -GF GT \
#'    -O output.table
#' }
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom readr read_csv
#' @importFrom dplyr filter tibble left_join
#'
#' @examples
#' \dontrun{
#' gatk_table_path <- "path/to/gatk_table.table"
#' col_data_path <- "path/to/col_data.csv"
#' drop_samples <- c("sample1", "sample2")
#' metadata <- list(population_structure = 'RIL', population_1_n = 20)
#'
#' bsa_experiment <- makeBSAExperimentFromGatkTable(gatk_table_path, col_data_path, drop_samples, metadata)
#' }
#'
#' @export
makeBSAExperimentFromGatkTable <- function(gatk_table_path,
                                           col_data_path = NULL,
                                           drop_samples = c(),
                                           metadata = list()) {
  # Function implementation
}

makeBSAExperimentFromGatkTable <- function(gatk_table_path,
                                           col_data_path,
                                           drop_samples = c(),
                                           metadata = list()) {

  .validate_makeBSAExperimentFromGatkTableInput(gatk_table_path,
                                                col_data_path,
                                                drop_samples,
                                                metadata)

  # read in the GATK table data. Note that this includes validation
  # of the necessary columns used below
  parsed_data = .read_in_gatk_table(gatk_table_path,
                                    drop_samples)

  # read in the col data, if it is provided
  if (!is.null(col_data_path)) {
    user_col_data = readr::read_csv(col_data_path, show_col_types = FALSE)
    # check if `samples` is a column in the col data
    if(!"sample" %in% colnames(user_col_data)){
      stop("The column `sample` must be in the col data.")
    }
    # check if all of the parsed_data$samples are in user_col_data$sample
    if(!all(parsed_data$samples %in% user_col_data$sample)){
      stop(sprintf("All of the samples in the GATK table must be in the col data. ",
                   "The following are missing: %s",
                   paste(base::setdiff(parsed_data$samples,
                                 user_col_data$sample),
                         collapse=", ")))
    }

    # Ensure that the set of samples in the user_col_data and the variant data
    # are the same
    user_col_data = user_col_data %>%
      dplyr::filter(sample %in% parsed_data$samples)
  }

  # create assay data list
  assays = .create_assay_list(parsed_data)

  # create rowRanges from the coordinate and variant level information
  # note that the rownames on assays are set in the .create_assays_sample
  # function
  row_ranges <- GenomicRanges::GRanges(seqnames = parsed_data$table$CHROM,
                                       ranges = IRanges::IRanges(start = parsed_data$table$POS,
                                                                 width = 1),
                                       strand = "*",
                                       ref = parsed_data$table$REF,
                                       alt = parsed_data$table$ALT,
                                       multiallelic = parsed_data$table$`MULTI-ALLELIC`,
                                       type = parsed_data$table$TYPE)

  col_data <- dplyr::tibble(sample = parsed_data$samples) %>%
    dplyr::left_join(user_col_data, by = "sample") %>%
    S4Vectors::DataFrame()
  rownames(col_data) = col_data$Sample

  # Create a new BSAExperiment object with the extracted data
  BSAExperiment(assays = assays,
                rowRanges = row_ranges,
                colData = col_data,
                metadata = metadata)
}


#' @importFrom stringr str_ends
#' @keywords internal
.validate_makeBSAExperimentFromGatkTableInput = function(gatk_table_path,
                                                         col_data_path,
                                                         drop_samples,
                                                         metadata) {

  error_messages <- c() # Initialize an empty vector to store error messages

  if(!is.null(drop_samples) & !is.character(drop_samples)){
    error_messages <- c(error_messages,
                        "`drop_samples` must be a character vector or NULL")
  }

  if(!file.exists(gatk_table_path)) {
    error_messages <- c(error_messages,
                        sprintf("`gatk_table_path` %s does not exist", gatk_table_path))
  }

  # validate col_data_path ----
  if(!is.null(col_data_path) & !file.exists(col_data_path)) {
    error_messages <- c(error_messages,
                        sprintf("`col_data_path` %s does not exist", col_data_path))
  }

  if(!is.null(col_data_path) &&
     !(stringr::str_ends(tolower(col_data_path), ".csv") | stringr::str_ends(tolower(col_data_path), ".csv.gz"))) {
    error_messages <-  c(error_messages,
                        paste0("`col_data_path` must be a csv file with ",
                                sprintf("`.csv` as an extension. Verify that %s is a ", col_data_path),
                                "csv and change the extension."))
  }

  # validate that the metadata is a list ----
  if(!is.list(metadata)){
    error_messages <- c(error_messages,
                        "`metadata` must be a list")
  }

  # If there are any error messages, stop and show all errors
  if(length(error_messages) > 0) {
    stop(paste(error_messages, collapse = "\n"))
  }
}

#' Read in GATK Table
#'
#' This function reads a GATK table from a file and extracts the sample names.
#' It checks that the table contains the required columns and sample columns.
#'
#' @param gatk_table_path The path to the GATK table file.
#' @param keep_multiallelic A boolean indicating whether to keep multi-allelic
#'   variants. Note: multi_allelic locations are not currently supported, so
#'   don't bother setting this to TRUE for the time being. Default is FALSE.
#' @param drop_samples A character vector of sample names to drop from the GATK
#'   table. Default is an empty character vector.
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect regex
#'
#' @return A list containing the GATK table and the sample names.
#'
#' @examples
#' gatk_table_path <- system.file("extdata", "your_gatk_table.txt", package = "your_package_name")
#' result <- .read_in_gatk_table(gatk_table_path)
#' names(result)
#'
#' @keywords internal
.read_in_gatk_table = function(gatk_table_path,
                               drop_samples = c(),
                               keep_multiallelic = FALSE){
  if(keep_multiallelic == TRUE){
    error_messages <- c(error_messages,
                        "`keep_multiallelic` is not yet implemented. Please set to FALSE.")
  }

  # Read the GATK table into a data frame
  gatk_table <- read.table(gatk_table_path,
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           check.names=FALSE)

  # ensure that the following columns are present in the GATK table:
  # CHROM, POS, REF, ALT, MULTI-ALLELIC TYPE
  if(!all(c("CHROM", "POS", "REF", "ALT", "MULTI-ALLELIC", "TYPE", "QUAL") %in% colnames(gatk_table))){
    stop(paste0("The GATK table does not contain the required columns: ",
                "CHROM, POS, REF, ALT, MULTI-ALLELIC, TYPE, QUAL. ",
                "See ?makeBSAExperimentFromgatk_table for more information."))
  }

  # confirm that at least the following columns are present in the sample
  # columns; AD, DP, PL, GQ. we can just grep for these specifically,
  # if they exist at all, we can assume that the sample columns are present
  if(!any(grepl("AD", colnames(gatk_table))) |
     !any(grepl("DP", colnames(gatk_table))) |
     !any(grepl("PL", colnames(gatk_table))) |
     !any(grepl("GQ", colnames(gatk_table)))){
    stop(paste0("The GATK table does not contain the required sample columns: ",
                "GT, AD, DP, PL, GQ. ",
                "See ?makeBSAExperimentFromgatk_table for more information."))
  }

  # Extract the sample names from column names that contain a '.'. All sample
  # columns from GATK varianttotable will have a period. non sample columns
  # won't
  sample_names <- unique(sub("\\..*$", "",
                             colnames(gatk_table)[grepl("\\.", colnames(gatk_table))]))

  # raise an error if the sample names are null
  if(length(sample_names) < 1){
    stop("No sample names found in the GATK table")
  }

  # if drop samples is not null, then check if the sampels listed are in
  # the sample_names. If not, raise a warning, but continue. If so, remove
  # those samples from the sample_names and also from the gatk_table
  if(!is.null(drop_samples)){
    if(!all(drop_samples %in% sample_names)){
      warning(paste0("The following samples are not in the GATK table: ",
                     paste0(drop_samples[!drop_samples %in% sample_names],
                            collapse=", ")))
    } else {
      sample_names <- sample_names[!sample_names %in% drop_samples]
      gatk_table <- gatk_table[,!grepl(paste0("^(", paste0(drop_samples, collapse="|"), ")\\."),
                                       colnames(gatk_table))]
    }
  }

  if(!keep_multiallelic){
    gatk_table = gatk_table %>%
      dplyr::filter(stringr::str_detect(`MULTI-ALLELIC`,
                                        stringr::regex('false', ignore_case=TRUE)))
  }

  list(table=gatk_table, samples=sample_names)
}

#'
#' @importFrom purrr map
#' @importFrom S4Vectors SimpleList
#'
#' @keywords internal
.create_assay_list = function(parsed_data){
  # create a SimpleList object to store the assays
  assays <- S4Vectors::SimpleList()
  # the AD value is a comma separated string where the first value is the REF
  # depth and the second value is the ALT depth
  assays$AD <- .parse_ad_table(as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".AD")]))
  assays$DP <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".DP")])
  assays$PL <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".PL")])
  assays$GQ <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".GQ")])


  # set the colnames of all the assays to the sample names
  assays <- purrr::map(assays, function(x) {colnames(x) <- parsed_data$samples; return(x)})
  # set rownames of all the assays
  assay_rownames <- paste0(parsed_data$table$CHROM, ":",
                           parsed_data$table$POS, "-",
                           parsed_data$table$REF, "-",
                           parsed_data$table$ALT)

  # add the rownames tot he assays and return the result
  purrr::map(assays, function(x) {rownames(x) <- assay_rownames; return(x)})
}

#' Parse AD Table
#'
#' This function takes a matrix of AD data, splits each element at the comma,
#' extracts the second value, converts these values to numeric, and reshapes
#' the result into a matrix with the same number of rows as the original AD data.
#'
#' @param ad_data A matrix of AD data.
#'
#' @importFrom dplyr mutate_all
#' @importFrom stringr str_split_fixed str_detect
#'
#' @return A numeric matrix where each element is the second value from the
#' corresponding element in the original AD data.
#'
#' @examples
#' ad_data <- matrix(c("1,149", "240,3", "225,65", "236,56", "406,141"), nrow = 5)
#' .parse_ad_table(ad_data)
#'
#' @keywords internal
.parse_ad_table = function(ad_data){
  # Define a function to extract the second part of a string
  extract_second_part <- function(x) {
    parts <- stringr::str_split_fixed(x, ",", 2)
    if (any(stringr::str_detect(parts[,2],","))) {
      problematic_data_subset = ad_data[which(stringr::str_detect(parts[,2],","))[1:1],]
      stop(paste0("Error with data:\n```\n",
                  paste(capture.output(dput(problematic_data_subset)), collapse = "\n"),
                  "\n```\nPlease provide your data to the code maintainer, preferrably via github, for help debugging."))
    }
    as.numeric(parts[,2])
  }

  # Apply the function to each column of the table
  tryCatch(
    apply(ad_data, 2, extract_second_part),
    warning = function(w) {
      list(
        message = paste("Parsing the AD table generated warnings. ",
                        "You should review these. ",
                        sprintf("%s", w$message),
                        collapse="\n"),
        warning = w
      )
    },
    error = function(e) {
      stop(
        paste(
          "Error parsing AD table.",
          "Please ensure that the AD table is in the correct format.",
          sprintf("Original error:\n```\n%s\n```", e$message),
          collapse = "\n"
        ),
        call. = FALSE
      )
    }
  )
}

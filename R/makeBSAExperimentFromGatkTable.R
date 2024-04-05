#'
#' Create a BSAExperiment object from a GATK `varianttotable` Table
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom readr read_csv
#' @importFrom dplyr filter tibble left_join
#'
#' @export
makeBSAExperimentFromGatkTable <- function(gatk_table_path,
                                           col_data_path,
                                           drop_samples = c(),
                                           high_confidence_depth=10,
                                           high_confidence_alt_percentage=.9,
                                           keep_multiallelic=FALSE,
                                           high_confidence_pl=NULL,
                                           high_confidence_gq=NULL) {

  .validate_makeBSAExperimentFromGatkTableInput(gatk_table_path,
                                                col_data_path,
                                                drop_samples,
                                                high_confidence_depth,
                                                high_confidence_alt_percentage,
                                                keep_multiallelic,
                                                high_confidence_pl,
                                                high_confidence_gq)

  # read in the GATK table data. Note that this includes validation
  # of the necessary columns used below
  parsed_data = .read_in_gatk_table(gatk_table_path,
                                    keep_multiallelic,
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
  assays = .create_assay_list(parsed_data,
                              high_confidence_depth,
                              high_confidence_alt_percentage,
                              high_confidence_pl,
                              high_confidence_gq)

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
  BSAExperiment(assays = assays, rowRanges = row_ranges, colData = col_data)
}


#' @keywords internal
.validate_makeBSAExperimentFromGatkTableInput = function(gatk_table_path,
                                                         col_data_path,
                                                         drop_samples,
                                                         high_confidence_depth,
                                                         high_confidence_alt_percentage,
                                                         keep_multiallelic,
                                                         high_confidence_pl,
                                                         high_confidence_gq) {

  error_messages <- c() # Initialize an empty vector to store error messages

  # validate .read_in_data arguments ----
  if(!is.logical(keep_multiallelic)){
    error_messages <- c(error_messages,
                        "`keep_multiallelic` must be a boolean")
  }

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
  if(!is.null(col_data_path) && tools::file_ext(col_data_path) != 'csv'){
    error_messages <- c(error_messages,
                        paste0("`col_data_path` must be a csv file with ",
                                sprintf("`.csv` as an extension. Verify that %s is a ", col_data_path),
                                "csv and change the extension."))
  }

  # validate assay related arguments ----
  if(!is.numeric(high_confidence_depth) | high_confidence_depth < 0) {
    error_messages <- c(error_messages,
                        "`high_confidence_depth` must be a positive integer")
  }
  if(!is.numeric(high_confidence_alt_percentage) | high_confidence_alt_percentage < 0 | high_confidence_alt_percentage > 1) {
    error_messages <- c(error_messages,
                        "`high_confidence_alt_percentage` must be a number between 0 and 1")
  }
  if(!is.null(high_confidence_pl) && (!is.numeric(high_confidence_pl) | high_confidence_pl < 0)) {
    error_messages <- c(error_messages,
                        "`high_confidence_pl` must be a positive number or NULL")
  }
  if(!is.null(high_confidence_gq) && (!is.numeric(high_confidence_gq) | high_confidence_gq < 0)) {
    error_messages <- c(error_messages,
                        "`high_confidence_gq` must be a positive number or NULL")
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
#' @param keep_multi_allelic A boolean indicating whether to keep multi-allelic
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
                               keep_multiallelic = FALSE,
                               drop_samples = c()){

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
  if(!all(c("CHROM", "POS", "REF", "ALT", "MULTI-ALLELIC", "TYPE") %in% colnames(gatk_table))){
    stop(paste0("The GATK table does not contain the required columns: ",
                "CHROM, POS, REF, ALT, MULTI-ALLELIC, TYPE. ",
                "See ?makeBSAExperimentFromgatk_table for more information."))
  }

  # confirm that at least the following columns are present in the sample
  # columns; GT, AD, DP, PL, GQ. we can just grep for these specifically,
  # if they exist at all, we can assume that the sample columns are present
  if(!any(grepl("GT", colnames(gatk_table))) |
     !any(grepl("AD", colnames(gatk_table))) |
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
.create_assay_list = function(parsed_data,
                              high_confidence_depth,
                              high_confidence_alt_percentage,
                              high_confidence_pl,
                              high_confidence_gq){
  # create a SimpleList object to store the assays
  assays <- S4Vectors::SimpleList()
  # the AD value is a comma separated string where the first value is the REF
  # depth and the second value is the ALT depth
  assays$GT <- .parse_gt_table(as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".GT")]),
                               parsed_data$table$REF)
  assays$AD <- .parse_ad_table(as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".AD")]))
  assays$DP <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".DP")])
  assays$PL <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".PL")])
  assays$GQ <- as.matrix(parsed_data$table[, paste0(parsed_data$samples, ".GQ")])
  assays$alt_percentage <- assays$AD / assays$DP
  # create a boolean matrix where the value is 1 if the depth is greater than
  # 10 and the alt_percentage is greater than .9 or less than .1
  # if high_confidence_pl and/or high_confidence_gq are provided, use those
  # to label high confidence calls, also
  assays$high_confidence <- (assays$DP >= high_confidence_depth) &
    ((assays$alt_percentage >= high_confidence_alt_percentage) |
       (assays$alt_percentage <= (1 - high_confidence_alt_percentage)))
  if(is.numeric(high_confidence_pl)){
    assays$high_confidence <- assays$high_confidence & (assays$PL <= high_confidence_pl)
  }
  if(is.numeric(high_confidence_gq)){
    assays$high_confidence <- assays$high_confidence & (assays$GQ >= high_confidence_gq)
  }
  # replace any NA values in high_confidence with FALSE
  assays$high_confidence[is.na(assays$high_confidence)] <- FALSE

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

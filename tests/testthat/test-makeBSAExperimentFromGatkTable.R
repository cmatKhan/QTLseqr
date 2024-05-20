test_that("test .validate_makeBSAExperimentFromGatkTable valid input", {
  mockery::stub(.validate_makeBSAExperimentFromGatkTableInput, "file.exists", TRUE)

  expect_silent(.validate_makeBSAExperimentFromGatkTableInput(
    gatk_table_path = "path/to/nonexistent_gatk_table.txt",
    col_data_path = "path/to/nonexistent_col_data.csv",
    drop_samples = NULL,
    metadata = list()
  ))
})

test_that("test .validate_makeBSAExperimentFromGatkTable error conditions", {
  mockery::stub(.validate_makeBSAExperimentFromGatkTableInput, "file.exists", FALSE)

  expect_error(
    .validate_makeBSAExperimentFromGatkTableInput(
      gatk_table_path = "nonexistent_gatk_table.txt", # Invalid because it does not exist
      col_data_path = "nonexistent_col_data.txt", # Invalid because it does not exist and wrong extension
      drop_samples = 123, # Invalid because it's not a character vector or NULL
      metadata = list()
    ),
    # Check for a comprehensive error message containing all validation issues
    paste(
      "`drop_samples` must be a character vector or NULL",
      "`gatk_table_path` nonexistent_gatk_table.txt does not exist",
      "`col_data_path` nonexistent_col_data.txt does not exist" ,
      "`col_data_path` must be a csv file with `.csv` as an extension. Verify that nonexistent_col_data.txt is a csv and change the extension.",
      sep = "\n"
    )
  )
})


test_that("test .read_in_gatk_table function", {

  file_paths = bsa_file_fixtures()
  gatk_table_path = file_paths$gatk_table
  # Call the function
  result <- .read_in_gatk_table(gatk_table_path)
  # Check that the result is a list with two elements: 'table' and 'samples'
  expect_type(result, "list")
  expect_equal(names(result), c("table", "samples"))

  # Check that 'table' is a data frame with the correct columns
  expect_s3_class(result$table, "data.frame")
  expect_equal(colnames(result$table),
               c("CHROM", "POS", "REF", "ALT", "MULTI-ALLELIC", "TYPE", "QUAL",
                 paste("C8", c('AD', 'DP', 'PL', 'GQ', 'GT'), sep='.'),
                 paste("KN99a", c('AD', 'DP', 'PL', 'GQ', 'GT'), sep='.'),
                 paste("SLB0021", c('AD', 'DP', 'PL', 'GQ', 'GT'), sep='.'),
                 paste("SLB0025", c('AD', 'DP', 'PL', 'GQ', 'GT'), sep='.')))

  # Check that 'samples' is a character vector with the correct values
  expect_type(result$samples, "character")
  expect_equal(result$samples, c("C8", "KN99a", 'SLB0021', 'SLB0025'))

  # Test the function with the 'drop_samples' parameter
  result_with_drop <- .read_in_gatk_table(gatk_table_path, drop_samples = c("C8"))
  expect_equal(result_with_drop$samples, c("KN99a", "SLB0021", "SLB0025"))
  expect_false("C8.AD" %in% colnames(result_with_drop$table))
})

test_that("parse_ad_table returns correct output", {

  ad_data = matrix(c("1,149","3,159","1,147","10,95","0,23",
                     "240,3", "237,0","223,0","42,0","28,0"),
                   nrow=5)
  colnames(ad_data) = c('sample1', 'sample2')

  result = .parse_ad_table(ad_data)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(5, 2))

  expected_result =   ad_data = matrix(c(149,159,147,95,23,
                                         3,0,0,0,0),
                                       nrow=5)
  colnames(expected_result) = c('sample1', 'sample2')

  expect_equal(result, expected_result)
})

test_that("test makeBSAExperimentFromGatkTable constructs an appropriate BSAExperiment object", {
  file_paths = bsa_file_fixtures()
  gatk_table_path = file_paths$gatk_table
  coldata_path = file_paths$coldata

  # Expected setup
  expected_samples <- c("SLB0021", "SLB0025")
  expected_chromosomes <- c("CP022321.1")
  expected_positions <- c(222, 558, 1443, 10528, 12097)

  file_paths = bsa_file_fixtures()
  gatk_table_path = file_paths$gatk_table
  coldata_path = file_paths$coldata

  expected_samples <- c("SLB0021", "SLB0025")
  expected_chromosomes <- c("CP022321.1")
  expected_positions <- c(222, 558, 1443, 10528, 12097)


  # Perform the operation and store result in the container
  result <- makeBSAExperimentFromGatkTable(
    gatk_table_path = gatk_table_path,
    col_data_path = coldata_path,
    drop_samples = c('KN99a', 'C8'),
    metadata = list(
      population_1_n = 200,
      population_2_n = 100,
      population_structure = 'RIL'))

  # Test if the result is a BSAExperiment object
  expect_true(inherits(result, "BSAExperiment"))

  # Validate rowRanges
  row_ranges <- rowRanges(result)
  expect_equal(as.character(seqnames(row_ranges)@values), 'CP022321.1')
  expect_equal(start(row_ranges), expected_positions)

  # Validate samples in colData
  col_data <- colData(result)
  expect_equal(col_data$sample, expected_samples)

  # Validate metadata
  expect_equal(result@metadata$population_1_n, 200)
  expect_equal(result@metadata$population_2_n, 100)
  expect_equal(result@metadata$population_structure, 'RIL')

})

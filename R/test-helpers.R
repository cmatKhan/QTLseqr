bsa_file_fixtures <- function(env = parent.frame()){

  gatk_table_path = withr::local_tempfile(fileext=".txt",
                                          .local_envir = env)
  coldata_path = withr::local_tempfile(fileext='.csv',
                                       .local_envir = env)

  # Create the GATK table in code
  gatk_table <- data.frame(
    CHROM = rep("CP022321.1", 5),
    POS = c(222, 558, 1443, 10528, 12097),
    REF = c("G", "C", "C", "A", "T"),
    ALT = c("A", "G", "T", "T", "A"),
    `MULTI-ALLELIC` = rep(FALSE, 5),
    TYPE = rep("SNP", 5),
    C8.AD = c("1,149", "3,159", "1,147", "10,95", "0,23"),
    C8.DP = c(150, 162, 148, 105, 23),
    C8.PL = c("4205,0", "4977,0", "4533,0", "1778,0", "490,0"),
    C8.GQ = c(99, 99, 99, 99, 99),
    C8.GT = c("A", "G", "T", "T", "A"),
    KN99a.AD = c("240,3", "237,0", "223,0", "42,0", "28,0"),
    KN99a.DP = c(243, 237, 223, 42, 28),
    KN99a.PL = c("0,7042", "0,7202", "0,6743", "0,1176", "0,700"),
    KN99a.GQ = c(99, 10, 6, 99, 99),
    KN99a.GT = c("G", "C", "C", "A", "T"),
    SLB0021.AD = c("225,65", "243,70", "230,61", "40,19", "31,8"),
    SLB0021.DP = c(290, 313, 292, 59, 39),
    SLB0021.PL = c("0,6352", "0,6879", "0,6604", "0,874", "0,620"),
    SLB0021.GQ = c(99, 10, 6, 99, 99),
    SLB0021.GT = c("G", "C", "C", "A", "T"),
    SLB0025.AD = c("84,17", "153,35", "143,24", "27,10", "54,4"),
    SLB0025.DP = c(101, 188, 167, 37, 58),
    SLB0025.PL = c("0,2666", "0,4700", "0,4721", "0,677", "0,1396"),
    SLB0025.GQ = c(99, 10, 6, 99, 99),
    SLB0025.GT = c("G", "C", "C", "A", "T"),
    check.names = FALSE
  )
  # Save the GATK table to a temporary file
  write.table(gatk_table, gatk_table_path, sep = "\t", row.names = FALSE)

  meta_df = data.frame(
    sample = c("SLB0021", "SLB0025"),
    group = c("ir7all", "ir7all"),
    pool = c(1, 1),
    replicate = c(1, 1),
    condition = c("inoculum", "Lung"),
    sac_day = c(0, 9),
    culture_time = c(0, 0)
  )

  write.table(meta_df, coldata_path, sep = ",", row.names = FALSE)

  return(list(gatk_table=gatk_table_path, coldata=coldata_path))
}

bsae_obj_fixture = function(add_comparisons = FALSE){

  paths = bsa_file_fixtures()

  bsae = makeBSAExperimentFromGatkTable(
    gatk_table_path = paths$gatk_table,
    col_data_path = paths$coldata,
    drop_samples = c('KN99a', 'C8'),
    high_confidence_depth = 10,
    high_confidence_alt_percentage = 0.9,
    keep_multiallelic = FALSE,
    high_confidence_pl = NULL,
    high_confidence_gq = NULL
  )

  if(add_comparisons){
    bsae = createComparisonsFrame(bsae,
                                  grouping_variable = "group",
                                  population_variable = "condition",
                                  population_base_condition = "inoculum",
                                  var1_name = 'population_1',
                                  var2_name = 'population_2',
                                  base_cond_in_each_group = TRUE)
  }

  bsae

}

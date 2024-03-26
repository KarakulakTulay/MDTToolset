library(MDTToolset)
library(data.table)
setup_cancer <- function() {
  data_matrix <- data.frame(
    ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004", "ENST000005", "ENST000006", "ENST000007"),
    ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002", "ENSG000003", "ENSG000003", "ENSG000003"),
    Sample1 = c(100, 200, 20, 50, 0, 40, 100),
    Sample2 = c(150, 250, 5, 10, 3, 20, 60),
    Sample3 = c(150, 250, 5, 10, 3, 20, 60),
    Sample4 = c(150, 250, 5, 10, 3, 20, 80),
    Sample5 = c(150, 250, 5, 10, 3, 20, 100),
    stringsAsFactors = FALSE
  )
  colnames(data_matrix)[3:7] <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  return(data_matrix)
}

setup_gtex <- function() {
  data_matrix <- data.frame(
    ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004", "ENST000005", "ENST000006", "ENST000007"),
    ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002", "ENSG000003", "ENSG000003", "ENSG000003"),
    Sample1 = c(100, 200, 50, 20, 0, 40, 20),
    Sample2 = c(150, 250, 5, 10, 5, 20, 10),
    Sample3 = c(150, 250, 5, 10, 5, 20, 10),
    Sample4 = c(150, 250, 5, 10, 5, 20, 8),
    Sample5 = c(150, 250, 5, 10, 5, 20, 7),
    stringsAsFactors = FALSE
  )
  colnames(data_matrix)[3:7] <-  c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  return(data_matrix)
}

test_that("switch_calculator returns correct output for simple input", {
  # Setup mock input data
  exp_values_pcawg <- setup_cancer()
  exp_values_gtex <- setup_gtex()
  cancer_ori <- find_mdts(exp_values_pcawg, 2, 2)
  gtex_ori <- find_mdts(exp_values_gtex, 2, 2)
  ensp_sequences <- data.frame(
    ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004", "ENST000005", "ENST000006", "ENST000007"),
    ENST_Seq = c("KLMAMGADEVRK", "KLMAMGADEVRALSP", "KLMAMGADEVRALAP", "KLMAMGADEVRALSKL",  "KLMAMGADEVRALAKL", "KLMAMGADEVRALKLKL",  "KLMAMGAAEVRALSKL")
  )
  cutoff <- 0
  cutoff_rate <- 2
  cutoff_other_MDTs <- 50
  ensp_only <- data.frame(
    Transcript.stable.ID = c("ENST000001", "ENST000002", "ENST000003", "ENST000004", "ENST000005", "ENST000006", "ENST000007")
  )

  # Expected output setup
  expected_output <- data.table(
    SampleID = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5"),
    ENSG = c("ENSG000003", "ENSG000003", "ENSG000003", "ENSG000003",  "ENSG000003"),
    dMDT = c("ENST000007", "ENST000007", "ENST000007", "ENST000007", "ENST000007"),
    ENST2_cancer = rep("ENST000006", 5),
    TPM1_cancer = c(100, 60, 60, 80, 100),
    TPM2_cancer = c(40, 20, 20, 20, 20),
    enrichment = c(2.5, 3.0, 3.0, 4.0, 5.0),
    p_value = rep(0.03125, 5),
    relative_cancer_exp = c(0.7142857, 0.7228916, 0.7228916, 0.7766990, 0.8130081),
    relative_gtex_exp = rep(0.2857143, 5),
    MDT_GTEx = rep("ENST000006", 5),
    adj_p = rep(0.03125, 5))

  # Call switch_calculator
  result <- switch_calculator(cancer_ori, gtex_ori, exp_values_pcawg, exp_values_gtex, ensp_sequences, cutoff, cutoff_rate, cutoff_other_MDTs, ensp_only)

  # Verify output
  expect_equal(result, expected_output, tolerance = 1e-7)
})

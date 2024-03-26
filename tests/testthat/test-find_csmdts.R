library(MDTToolset)
library(data.table)

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

exp_values_gtex <- setup_gtex()
sample_mdt_file <- find_mdts(exp_values_gtex, 2, 2)

test_that("find_csmdts filters data correctly based on cutoff", {
  # Test with a cutoff that filters out some data
  result <- find_csmdts(sample_mdt_file, 50) # Assuming 50% as cutoff
  expected_ENST1 <- c("ENST000004", 'ENST000006')
  expect_equal(sort(unique(result$ENST1)), sort(expected_ENST1))
  # Test to ensure the function returns dataframe
  expect_s3_class(result, "data.frame")
})

test_that("find_csmdts filters data correctly based on cutoff 100", {
  # Test with a cutoff that filters out some data
  result <- find_csmdts(sample_mdt_file, 100) # Assuming 50% as cutoff
  expected_ENST1 <- c("ENST000006")
  expect_equal(sort(unique(result$ENST1)), sort(expected_ENST1))
  # Test to ensure the function returns dataframe
  expect_s3_class(result, "data.frame")
})

test_that("find_csmdts handles empty input correctly", {
  # Test with empty data frame
  empty_data <- data.frame(SampleID = character(), ENST1 = character(), stringsAsFactors = FALSE)
  result <- find_csmdts(empty_data, 50)
  expect_equal(nrow(result), 0)
  expect_s3_class(result, "data.frame")
})






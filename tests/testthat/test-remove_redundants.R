library(MDTToolset)
library(data.table)

setup <- function() {
  data_matrix <- data.frame(
    ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004"),
    ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002"),
    Sample1 = c(100, 100, 20, 50),
    Sample2 = c(100, 100, 5, 10),
    Sample3 = c(100, 100, 5, 10),
    Sample4 = c(100, 100, 5, 10),
    Sample5 = c(100, 100, 5, 10),
    stringsAsFactors = FALSE
  )
  colnames(data_matrix)[3:7] <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  return(data_matrix)
}

# Test preparation: Load or create the necessary objects for testing
data_matrix <- setup() # Create a sample data matrix
enst_ensp_sequences <- data.frame(
  ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004"),
  ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002"),
  ENST_Seq = c('GGGACAGGGGGC', 'GGGACAGGGGGC', 'GAAATAGT', 'GAAATAGTS'),
  stringsAsFactors = FALSE
) # Sample mane selects
mane_select <- data.frame(
  ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002"),
  ENSG.s = c("ENSG000001.1", "ENSG000001.1", "ENSG000002.1", "ENSG000002.1"),
  ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004"),
  ENST.s = c("ENST000001.1", "ENST000002.1", "ENST000003.1", "ENST000004.1"),
  Mane = c('Mane', 'WithoutMane', 'Mane', 'WithoutMane')
) # Create a sample mane transcript list

# Test 1: Verify that function returns a dataframe
test_that("remove_redundants returns a dataframe", {
  result <- remove_redundants(data_matrix, enst_ensp_sequences, mane_select)
  expect_s3_class(result, "data.frame")
})

# Test 2: Verify column names after processing
test_that("column names are correctly formatted", {
  result <- remove_redundants(data_matrix, enst_ensp_sequences, mane_select)
  expected_colnames <- c("ENST", "ENSG", "Sample1", "Sample2", "Sample3", "Sample4", "Sample5") # Specify the expected column names after processing
  expect_equal(colnames(result), expected_colnames)
})

# Test 3: Check for the removal of redundant transcripts
test_that("redundant transcripts are removed", {
  result <- remove_redundants(data_matrix, enst_ensp_sequences, mane_select)
  # Assuming you have a way to calculate the expected number of unique transcripts
  expected_unique_transcripts <- c('ENST000001', 'ENST000003', 'ENST000004')
  observed_unique_transcripts <- unique(result$ENST)
  expect_equal(observed_unique_transcripts, expected_unique_transcripts)
})

# Test 4: Ensure TPMs are multiplied correctly for non-redundant transcripts
test_that("TPMs are multiplied correctly", {
  specific_data_matrix <- setup()
  result <- remove_redundants(specific_data_matrix, enst_ensp_sequences, mane_select)

  # Verify the TPM values are as expected
  expected_tpm_values <- data_matrix[data_matrix$ENST == 'ENST000001', 3:7]*2 # Specify expected TPM values
  observed_tpm_values <- result[result$ENST == 'ENST000001', 3:7]
  expect_equal(observed_tpm_values, expected_tpm_values)
})



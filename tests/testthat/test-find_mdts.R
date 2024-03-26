library(MDTToolset)

# Setup mock data for testing
setup <- function() {
  data_matrix <- data.frame(
    ENST = c("ENST000001", "ENST000002", "ENST000003", "ENST000004"),
    ENSG = c("ENSG000001", "ENSG000001", "ENSG000002", "ENSG000002"),
    Sample1 = c(100, 200, 0, 50),
    Sample2 = c(150, 250, 10, 5),
    stringsAsFactors = FALSE
  )
  colnames(data_matrix)[3:4] <- c("Sample1_TPM", "Sample2_TPM")
  return(data_matrix)
}

# Test that it returns a dataframe
test_that("find_mdts returns a dataframe", {
  data_matrix <- setup()
  result <- MDTToolset::find_mdts(data_matrix, cutoff = 1.5, min_exp = 50)
  expect_true(is.data.frame(result))
})

# Test for correct number of columns in the output
test_that("Output has the correct number of columns", {
  data_matrix <- setup()
  result <- MDTToolset::find_mdts(data_matrix, cutoff = 1.5, min_exp = 50)
  expect_equal(ncol(result), 7)
})

# Test for correct filtering based on min_exp and cutoff
test_that("Function filters based on min_exp and cutoff correctly", {
  data_matrix <- setup()
  result <- MDTToolset::find_mdts(data_matrix, cutoff = 1.5, min_exp = 50)
  # Verify if all returned rows comply with min_exp and cutoff criteria
  expect_true(all(result$TPM1 >= 50 & result$rate >= 1.5))
})


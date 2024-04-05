#' Find MDTs in each sample
#'
#' @param data_matrix non-redundant expression data of transcripts
#' @param cutoff minimum enrichment value
#' @param min_exp minimum cutoff value for transcript expression
#'
#' @return dataframe
#' @importFrom dplyr group_by mutate select count n arrange filter distinct desc
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'

find_mdts <- function(data_matrix, cutoff, min_exp) {
  out2 <- NULL
  data_matrix <- data_matrix %>%
    dplyr::group_by(.data$ENSG) %>%
    dplyr::filter(dplyr::n() > 1)

  for (eachcol in 3:ncol(data_matrix)) {
    whole_data <- NULL
    colname <- colnames(data_matrix)[eachcol]
    new_data <- data_matrix %>%
      dplyr::group_by(.data$ENSG) %>%
      dplyr::arrange(desc(dplyr::across(all_of(colname))), .by_group = TRUE) %>%
      dplyr::filter(dplyr::row_number() %in% 1:2) %>%
      dplyr::select(1, 2, .data[[colname]]) %>%
      as.data.frame()

    ENSTs <- new_data[, c(1, 2)] %>%
      dplyr::group_by(.data$ENSG) %>%
      dplyr::mutate(group_index = 1:dplyr::n() %>% paste0("enst_", .)) %>%
      tidyr::pivot_wider(names_from = group_index, names_prefix = "ENST_", values_from=ENST,  names_sep = "_")

    TPMs <- new_data[, c(2, 3)] %>%
      dplyr::group_by(.data$ENSG) %>%
      dplyr::mutate(group_index = 1:dplyr::n() %>% paste0("tpm_", .)) %>%
      tidyr::pivot_wider(names_from = group_index, names_prefix = "TPM", values_from=colname,  names_sep = "_")

    whole_data <- cbind(colname, ENSTs, TPMs[, c(2, 3)])

    colnames(whole_data) <- c("SampleID", "ENSG", "ENST1", "ENST2", "TPM1", "TPM2")

    whole_data$rate <- ifelse(whole_data$TPM2 == 0 | whole_data$ENST2 == "NA", whole_data$TPM1, whole_data$TPM1 / whole_data$TPM2)

    out <- whole_data %>% dplyr::filter(.data$TPM1 >= min_exp & .data$rate >= cutoff)
    out2 <- rbind(out2, out)
  }

  return(as.data.frame(out2))
}

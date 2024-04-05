#' Remove Redundant Transcripts
#'
#' @param data_matrix Bulk RNA Seq expression matrix from samples, columns are samples and first row is gene ids, second row is transcript ids and the rest of the rows are samples
#' @param enst_ensp_sequences amino acid sequence of each ensp
#' @param mane_select mane transcript list is needed to decide which transcript will be kept
#'
#' @return dataframe
#' @importFrom dplyr group_by mutate select count n arrange filter distinct desc everything
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

remove_redundants <- function(data_matrix, enst_ensp_sequences, mane_select) {
  # remove redundant transcripts and multiply other transcript TPMs with the number of redundant transcripts
  change_dot <- function(x) {
    x <- as.character(x)
    x <- gsub(pattern = ".", replacement = "-", x = x, fixed = TRUE)
    return(x)
  }
  colnames(data_matrix) <- change_dot(colnames(data_matrix))
  colnames(data_matrix) <- gsub("^X", "", colnames(data_matrix))

  mane_select_essentials <- mane_select[, c(1, 2, 4)]
  colnames(mane_select_essentials) <- c("ENSG", "ENST", "Mane")

  enst_ensp_sequences_mane_select <- merge(enst_ensp_sequences, mane_select_essentials, all.x = TRUE, sort = FALSE)

  enst_ensp_sequences_mane_select$Mane <- ifelse(enst_ensp_sequences_mane_select$Mane != "NA", "Mane", "WithoutMane")

  cleaned_list <- enst_ensp_sequences_mane_select %>%
    dplyr::arrange(.data$ENST, desc(.data$Mane == "Mane")) %>%
    dplyr::distinct(.data$ENSG, .data$ENST_Seq, .keep_all = TRUE)

  data_matrix_cleaned <- data_matrix[data_matrix$ENST %in% cleaned_list$ENST, ]

  counts_ensgs <- enst_ensp_sequences %>%
    dplyr::group_by(.data$ENSG, .data$ENST_Seq) %>%
    dplyr::reframe(n = dplyr::n(), dplyr::across(.cols = everything(), ~ .)) %>%
    dplyr::arrange(desc(n))

  counts_ensgs_cleaned <- counts_ensgs[counts_ensgs$ENST %in% data_matrix_cleaned$ENST, ]

  ordered_data_matrix <- merge(data_matrix_cleaned, counts_ensgs_cleaned[, c(3,4)]) %>% dplyr::distinct()

  ordered_data_matrix_multipyled <- ordered_data_matrix[,3:(ncol(ordered_data_matrix)-1)]*ordered_data_matrix[,ncol(ordered_data_matrix)]

  ordered_data_matrix_multipyled_with_names <- cbind(ordered_data_matrix[, 1:2], ordered_data_matrix_multipyled)

  return(ordered_data_matrix_multipyled_with_names)
}

#' MDTs found in %n of samples
#'
#' @param samples_mdt_files mdt files that is found by find_mdt() function
#' @param cutoff cutoff to define percentages of samples having that MDT in the cohort
#'
#' @return dataframe
#' @importFrom dplyr group_by mutate select count n arrange filter distinct desc
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

find_csmdts <- function(samples_mdt_files, cutoff) {
  gtex_wo_redundant_tissue <- NULL
  numberofsamples <- NULL
  data_with_percentages <- NULL
  data_after_cutoff <- NULL

  numberofsamples <- length(unique(samples_mdt_files$SampleID))
  data_with_percentages <- samples_mdt_files %>%
    dplyr::group_by(.data$ENST1) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(freq = .data$n / numberofsamples * 100) %>%
    dplyr::arrange(desc(.data$freq))

  data_after_cutoff <- data_with_percentages %>% dplyr::filter(.data$freq >= cutoff)

  gtex_wo_redundant_tissue <- samples_mdt_files %>% dplyr::filter(.data$ENST1 %in% data_after_cutoff$ENST1)

  return(gtex_wo_redundant_tissue)
}

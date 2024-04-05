#' Plot MDT Distribitions
#'
#' @param cond1_MDTs  output of find_mdts function for the condition 1
#' @param cond2_MDTs  output of find_mdts function for the condition 2
#'
#' @return ggplot
#' @importFrom ggplot2 geom_boxplot
#' @importFrom dplyr group_by summarize across arrange distinct count
#' @export

MDT_boxplot <- function(cond1_MDTs, cond2_MDTs, title, x_axis_label, y_axis_label) {

  MDT_counts_cond1 <- cond1_mdts %>% dplyr::group_by(SampleID, ENSG) %>% dplyr::summarize(n=n(), across()) %>% arrange(desc(n)) %>% distinct() %>% dplyr::count(SampleID)
  MDT_counts_cond1$Type <- 'cond1'
  MDT_counts_cond2 <- cond2_mdts %>% dplyr::group_by(SampleID, ENSG) %>% dplyr::summarize(n=n(), across()) %>% arrange(desc(n)) %>% distinct() %>% dplyr::count(SampleID)
  MDT_counts_cond2$Type <- 'cond2'


  title = title
  x_axis_label = x_axis_label
  y_axis_label = y_axis_label
  MDT_counts <- rbind(MDT_counts_cond1, MDT_counts_cond2)


  mdt_boxplot <- ggplot2::ggplot(MDT_counts, aes(y=n, x = Type, fill=Type)) +

    geom_boxplot(outlier.shape=8,outlier.colour="red", fill='#f7fcf0') +

    labs(
      title = title,
      x = x_axis_label,
      y = y_axis_label
    ) +

    theme(
      plot.title = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14, angle=45),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14)) + theme_bw()

  return(mdt_boxplot)
}


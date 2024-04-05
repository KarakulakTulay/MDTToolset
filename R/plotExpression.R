#' Title
#'
#' @param dMDT ENST id
#' @param MDT ENST id
#' @param sampleName sample id
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise ungroup
#' @return plot
#' @export
#'
plotExpression <- function(dMDT, MDT, sampleName, cond1_exp, cond2_exp) {

  dMDT_TPM_cond1 <- cond1_exp[cond1_exp$ENST == dMDT ,sampleName]
  MDT_TPM_cond1 <- cond1_exp[cond1_exp$ENST == MDT,sampleName]
  dMDT_TPM_cond2 <- as.vector(unlist(cond2_exp[cond2_exp$ENST == dMDT, 3:ncol(cond2_exp)]))
  MDT_TPM_cond2 <- as.vector(unlist(cond2_exp[cond2_exp$ENST == MDT, 3:ncol(cond2_exp)]))

  dMDT_MDT_TPM <- data.frame(names = c('dMDT', 'MDT', rep('dMDT', length(dMDT_TPM_cond2)), rep('MDT', length(MDT_TPM_cond2))), TPM = c(dMDT_TPM_cond1, MDT_TPM_cond1, dMDT_TPM_cond2, MDT_TPM_cond2),
                             Condition = c('cond1', 'cond1',  rep('cond2', length(dMDT_TPM_cond2)),  rep('cond2', length(dMDT_TPM_cond2))))

  # Calculate mean and standard error for cond2
  cond2_stats <- dMDT_MDT_TPM %>%
    dplyr::filter(Condition == "cond2") %>%
    dplyr::group_by(names) %>%
    dplyr::summarise(mean = mean(TPM),
                     se = sd(TPM) / sqrt(n())) %>%
    dplyr::ungroup()

  colnames(cond2_stats) <- c('names', 'TPM', 'se')
  cond1_stats <- dMDT_MDT_TPM[dMDT_MDT_TPM$Condition == 'cond1',1:2]
  cond1_stats$se <- 0
  #cond1_stats$sd <- 0
  cond2_stats$Condition <- 'cond2'
  cond1_stats$Condition <- 'cond1'

  cond_stats <- rbind(cond1_stats, cond2_stats)
  # Start plotting
  plot <- ggplot(data = cond_stats, aes(x = names, y = TPM, fill = names)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    facet_grid(~Condition) +
    labs(x = "dMDT and MDT", y = "TPM", title = "TPM values of dMDT and MDT by Condition") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#D95F02", "#1B9E77"))  + # Set your own colors
    geom_errorbar(data = cond2_stats, aes(x = names, ymin = TPM - se, ymax = TPM + se),
                position = position_dodge(width = 0.7), width = 0.25, inherit.aes = FALSE)


  return(plot)
}

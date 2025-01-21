plot_strand_bias_cus <- function (strand_bias, colors = NA, sig_type = c("fdr", "p")) 
{
  type <- significant <- strand_1 <- strand_2 <- log2_ratio <- NULL
  sig_plot <- log2_ratio_no1pseudo <- NULL
  sig_type <- match.arg(sig_type)
  if (.is_na(colors)) {
    # colors <- COLORS6
    colors <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  var_names <- colnames(strand_bias)[3:4]
  colnames(strand_bias)[3:4] <- c("strand_1", "strand_2")
  strand_bias <- dplyr::mutate(strand_bias, log2_ratio = log2((strand_1 + 
                                                                 0.1)/(strand_2 + 0.1)), log2_ratio_no1pseudo = log2((strand_1)/(strand_2 + 
                                                                                                                                   0.1)))
  if (sig_type == "p") {
    strand_bias$sig_plot <- strand_bias$significant
  }
  else {
    strand_bias$sig_plot <- strand_bias$significant_fdr
  }
  max <- round(max(abs(strand_bias$log2_ratio)), digits = 1) + 
    0.1
  pos_stars <- abs(strand_bias$log2_ratio_no1pseudo)
  max_pos_star <- round(max(pos_stars[is.finite(pos_stars)]), 
                        digits = 1) + 0.1
  if (max < max_pos_star) {
    max <- max_pos_star
  }
  label2 <- log2(strand_bias$ratio)
  select <- which(is.finite(label2))
  label2[select] <- " "
  plot <- ggplot(strand_bias, aes(x = type, y = log2_ratio, 
                                  fill = type)) + geom_bar(colour = "black", stat = "identity", 
                                                           position = "identity") + scale_fill_manual(values = c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")) + 
    scale_y_continuous(limits = c(-max, max)) + geom_text(aes(x = type, 
                                                              y = log2_ratio_no1pseudo, label = sig_plot, vjust = ifelse(sign(log2_ratio_no1pseudo) > 
                                                                                                                           0, 0.5, 1)), size = 8, position = ggplot2::position_dodge(width = 1)) + 
    facet_grid(. ~ fct_rev(group)) + theme_bw() + theme(axis.ticks = element_blank(), 
                                               axis.text.x = element_blank(), legend.title = element_blank()) + 
    xlab("") + ylab(paste("log2(", var_names[1], "/", var_names[2], 
                          ")", sep = "")) + scale_x_discrete(breaks = NULL)
  return(plot)
}

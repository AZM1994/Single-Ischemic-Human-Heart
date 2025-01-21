plot_strand_cus <- function (strand_bias_df, mode = c("relative", "absolute"), 
          colors = NA) 
{
  type <- relative_contribution <- no_mutations <- y_vals <- NULL
  mode <- match.arg(mode)
  if (.is_na(colors)) {
    # colors <- COLORS6
    colors <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  if (mode == "relative") {
    strand_bias_df$y_vals <- strand_bias_df$relative_contribution
    y_lab <- "Relative contribution"
  }
  else if (mode == "absolute") {
    strand_bias_df$y_vals <- strand_bias_df$no_mutations
    y_lab <- "Total number of mutations"
  }
  withCallingHandlers({
    plot <- ggplot(strand_bias_df, aes(x = type, y = y_vals, 
                                       fill = type, alpha = strand)) + geom_bar(stat = "identity", 
                                                                                position = "dodge", colour = "black", cex = 0.5) + 
      scale_fill_manual(values = colors) + scale_alpha_discrete(range = c(1, 
                                                                          0.4)) + labs(y = y_lab, x = "") + facet_grid(. ~ 
                                                                                                                         rev(group)) + theme_bw() + scale_x_discrete(breaks = NULL)
  }, warning = function(w) {
    if (grepl("Using alpha for a discrete variable is not advised.", 
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
  return(plot)
}

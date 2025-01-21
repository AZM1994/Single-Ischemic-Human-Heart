plot_bootstrapped_contribution_conditional <- function (contri_boots, mode = c("absolute", "relative"), plot_type = c("jitter", 
                                                                        "barplot", "dotplot")) 
{
  mode <- match.arg(mode)
  plot_type <- match.arg(plot_type)
  sig <- contri <- lower <- upper <- percentage <- NULL
  if (mode == "relative") {
    contri_boots <- contri_boots/rowSums(contri_boots)
    ylab_text <- "Relative mutation contribution"
    jitter_height <- 0.02
  }
  else {
    ylab_text <- "Nr. contributed mutations"
    jitter_height <- 0.2
  }
  contri_tb <- contri_boots %>% as.data.frame() %>% tibble::rownames_to_column("exp") %>% 
    tidyr::gather(key = "sig", value = "contri", -exp) %>% 
    dplyr::mutate(sample = gsub("_[^_]+$", "", exp), sample = factor(sample, 
                                                                     levels = unique(sample)), sig = factor(sig, levels = unique(sig)))
  nr_sigs <- length(unique(contri_tb$sig))
  if (plot_type == "jitter") {
    fig <- ggplot(contri_tb, aes(x = sig, y = contri, color = sig)) + 
      geom_jitter(stat = "identity", height = jitter_height, 
                  size = 0.3) + scale_color_discrete(guide = "none") + 
      facet_grid(sample ~ .) + labs(y = ylab_text)
  }
  else if (plot_type == "barplot") {
    contri_tb2 <- contri_tb %>% dplyr::group_by(sample, 
                                                sig) %>% dplyr::summarise(mean = mean(contri), lower = quantile(contri, 
                                                                                                                0.025), upper = quantile(contri, 0.975)) %>% dplyr::ungroup()
    fig <- ggplot(contri_tb2, aes(x = sig, y = mean, fill = sig)) + 
      geom_bar(stat = "identity") + geom_errorbar(aes(ymin = lower, 
                                                      ymax = upper), width = 0.2) + scale_fill_discrete(guide = "none") + 
      facet_grid(sample ~ .) + labs(y = ylab_text)
  }
  else if (plot_type == "dotplot") {
    contri_tb3 <- contri_tb %>% dplyr::group_by(sample, 
                                                sig) %>% dplyr::summarise(mean = mean(contri[contri != 
                                                                                               0]), percentage = sum(contri != 0)/dplyr::n()) %>% 
      dplyr::ungroup() %>% dplyr::filter(!is.na(mean)) %>% 
      dplyr::mutate(sample = factor(sample, levels = rev(levels(sample))))
    max_dot_size <- dplyr::case_when(nr_sigs >= 40 ~ 5, 
                                     nr_sigs >= 30 ~ 7, nr_sigs >= 20 ~ 8, nr_sigs >= 
                                       10 ~ 10, TRUE ~ 15)
    fig <- ggplot(contri_tb3, aes(x = sig, y = sample)) + 
      geom_point(aes(color = percentage, size = mean)) + 
      scale_color_distiller(palette = "RdYlBu", limits = c(0, 
                                                           1)) + scale_size_continuous(range = c(1, max_dot_size)) + 
      labs(size = "mean contribution", colour = "percentage contribution", 
           y = "Sample")
  }
  fig <- fig + labs(x = "Signature") + theme_classic() + theme(axis.text.x = element_text(angle = 90, 
                                                                                          size = 10, hjust = 1, vjust = 0.5), text = element_text(size = 12), 
                                                               strip.text.y = element_text(size = 8))
  if (plot_type == "dotplot") {
    fig <- fig + theme(panel.grid.major = element_line(colour = "gray92"))
  }
  return(fig)
}

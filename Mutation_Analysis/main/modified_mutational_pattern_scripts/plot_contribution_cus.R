plot_contribution_customize <- function (contribution, signatures = NA, index = NA, coord_flip = FALSE, 
          mode = c("relative", "absolute"), palette = NA, order_list = NA) 
{
  mode <- match.arg(mode)
  if (!.is_na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }
  Sample <- Contribution <- Signature <- NULL
  if (mode == "absolute" & !.is_na(signatures)) {
    total_signatures <- colSums(signatures)
    abs_contribution <- contribution * total_signatures
  }
  tb <- contribution %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
    tidyr::pivot_longer(-Signature, names_to = "Sample", 
                        values_to = "Contribution") %>% dplyr::mutate(Sample = factor(Sample, 
                                                                                      levels = unique(Sample)), Signature = factor(Signature, 
                                                                                                                                   levels = unique(Signature)))
  if (mode == "absolute") {
    bar_geom <- geom_bar(stat = "identity", colour = "black")
    y_lab <- "Absolute contribution \n (no. mutations)"
  }
  else if (mode == "relative") {
    bar_geom <- geom_bar(position = "fill", stat = "identity", 
                         colour = "black")
    y_lab <- "Relative contribution"
  }
  present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% unique()
  plot <- ggplot(tb, aes(x=factor(Sample, level=order_list), y = Contribution, fill = Signature)) + 
    bar_geom + labs(x = "", y = y_lab) + scale_fill_discrete(breaks = present_sigs) + 
    theme_bw() + theme(panel.grid.minor.x = element_blank(), 
                       panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
                       panel.grid.major.y = element_blank())
  if (!.is_na(palette)) {
    plot <- plot + scale_fill_manual(name = "Signature", 
                                     values = palette)
  }
  if (coord_flip) {
    plot <- plot + coord_flip() + xlim(rev(levels(factor(tb$Sample))))
  }
  else {
    plot <- plot + xlim(levels(factor(tb$Sample)))
  }
  return(plot)
}

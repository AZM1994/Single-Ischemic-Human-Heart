plot_96_profile_abs <- function (mut_matrix, colors = NA, ymax = NA, condensed = FALSE) 
{
  freq <- full_context <- substitution <- context <- NULL
  # ymax <- apply(X = mut_matrix, MARGIN = 2, FUN = max)
  # ymax <- max(mut_matrix)
  if (.is_na(colors)) {
    # colors <- COLORS6
    # colors <- c("#2EBAED", "#000000", "#DE1C14", "#E98C7B", "#D4D2D2", "#ADCC54", 
    #             "#F0D0CE")
    colors <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  # norm_mut_matrix <- apply(mut_matrix, 2, function(x) x)
  norm_mut_matrix <- mut_matrix
  # norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context, 
                                                      "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
                                                                                                                 "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
                        values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
                                                                              levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  my_breaks <- function(x){seq(min(x), max(x), length.out = 5)}
  my_labels <- function(x){round(x, digits = 2)}
  plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, 
                                width = width)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = colors) + 
    facet_grid(sample ~ substitution) + 
    # facet_grid(sample ~ substitution, scales="free") + 
    ylab("Absolute contribution") + 
    # coord_cartesian(ylim = c(0, NA)) +
    # scale_y_continuous(breaks = seq(0, ymax, 20)) +
    # scale_y_continuous(breaks = my_breaks, labels = my_labels) +
    guides(fill = "none") + theme_bw() + theme(axis.title.y = element_text(size = 12, 
                                                                           vjust = 1), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
                                               axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5), 
                                               strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9), 
                                               panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
                                                                                                            "lines"))
  return(plot)
}

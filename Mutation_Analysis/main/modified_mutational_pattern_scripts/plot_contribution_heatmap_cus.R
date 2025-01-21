plot_contribution_heatmap_cus <- function (contribution, sig_order = NA, sample_order = NA, cluster_samples = TRUE, 
          cluster_sigs = FALSE, method = "complete", plot_values = FALSE) 
{
  Signature <- Sample <- Contribution <- x <- y <- xend <- yend <- NULL
  if (!inherits(contribution, "matrix")) {
    stop("contribution must be a matrix")
  }
  if (is.null(row.names(contribution))) {
    stop("contribution must have row.names (signature names)")
  }
  contribution <- t(contribution)
  # contribution_norm <- contribution/rowSums(contribution)
  contribution_norm <- contribution
  if (!.is_na(sample_order) & cluster_samples == TRUE) {
    stop("sample_order can only be provided when cluster_samples is FALSE", 
         call. = FALSE)
  }
  else if (!.is_na(sample_order)) {
    if (!inherits(sample_order, "character")) {
      stop("sample_order must be a character vector", 
           call. = FALSE)
    }
    if (length(sample_order) != nrow(contribution_norm)) {
      stop("sample_order must have the same length as the number\n          of samples in the explained matrix", 
           call. = FALSE)
    }
  }
  else if (cluster_samples == TRUE) {
    hc.sample <- hclust(dist(contribution_norm), method = method)
    sample_order <- rownames(contribution_norm)[hc.sample$order]
    dhc <- as.dendrogram(hc.sample)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram_rows <- ggplot(ggdendro::segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + scale_y_reverse(expand = c(0.2, 0)) + 
      ggdendro::theme_dendro()
  }
  else {
    sample_order <- rownames(contribution_norm)
  }
  if (!.is_na(sig_order) & cluster_sigs == TRUE) {
    stop("sig_order can only be provided when cluster_sigs is FALSE", 
         call. = FALSE)
  }
  else if (!.is_na(sig_order)) {
    if (!inherits(sig_order, "character")) {
      stop("sig_order must be a character vector", call. = FALSE)
    }
    if (length(sig_order) != ncol(contribution_norm)) {
      stop("sig_order must have the same length as the number\n           of signatures in the explained matrix", 
           call. = FALSE)
    }
  }
  else if (cluster_sigs == TRUE) {
    hc.sample2 <- contribution_norm %>% t() %>% dist() %>% 
      hclust(method = method)
    sig_order <- colnames(contribution_norm)[hc.sample2$order]
    dhc <- as.dendrogram(hc.sample2)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram_cols <- ggplot(ggdendro::segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      ggdendro::theme_dendro() + scale_y_continuous(expand = c(0.2, 
                                                               0))
  }
  else {
    sig_order <- colnames(contribution_norm)
  }
  contribution_norm.m <- contribution_norm %>% as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% tidyr::pivot_longer(-Sample, 
                                                                 names_to = "Signature", values_to = "Contribution") %>% 
    dplyr::mutate(Signature = factor(Signature, levels = sig_order), 
                  Sample = factor(Sample, levels = sample_order))
  heatmap <- ggplot(contribution_norm.m, aes(x = Signature, 
                                             y = Sample, fill = Contribution, order = Sample)) + 
    geom_raster() + scale_fill_distiller(palette = "YlGnBu", 
                                         direction = 1, name = "Relative \ncontribution", limits = c(0, 
                                                                                                     500)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                                                                                                                         hjust = 1, vjust = 0.5), panel.grid.major = element_blank(), 
                                                                                                                              panel.grid.minor = element_blank()) + labs(x = NULL, 
                                                                                                                                                                         y = NULL)
  if (plot_values) {
    heatmap <- heatmap + geom_text(aes(label = round(Contribution, 
                                                     2)), size = 3)
  }
  if (cluster_samples == TRUE & cluster_sigs == TRUE) {
    empty_fig <- ggplot() + theme_void()
    plot_final <- cowplot::plot_grid(empty_fig, dendrogram_cols, 
                                     dendrogram_rows, heatmap, align = "hv", axis = "tblr", 
                                     rel_widths = c(0.3, 1), rel_heights = c(0.3, 1))
  }
  else if (cluster_samples == TRUE & cluster_sigs == FALSE) {
    plot_final <- cowplot::plot_grid(dendrogram_rows, heatmap, 
                                     align = "h", rel_widths = c(0.3, 1))
  }
  else if (cluster_samples == FALSE & cluster_sigs == TRUE) {
    plot_final <- cowplot::plot_grid(dendrogram_cols, heatmap, 
                                     align = "v", rel_heights = c(0.3, 1)) + ylim(rev(levels(factor(contribution_norm.m$Sample))))
  }
  else {
    plot_final <- heatmap + ylim(rev(levels(factor(contribution_norm.m$Sample))))
  }
  return(plot_final)
}

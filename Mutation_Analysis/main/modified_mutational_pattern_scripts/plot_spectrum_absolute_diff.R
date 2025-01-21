plot_spectrum_absolute_diff <- function (type_occurrences, CT = FALSE, by = NA, indv_points = FALSE, 
                                    error_bars = c("95%_CI", "stdev", "SEM", "none"), colors = NA, 
                                    legend = TRUE, condensed = FALSE) 
{
  value <- nmuts <- sub_type <- variable <- error_pos <- stdev <- total_mutations <- NULL
  x <- total_individuals <- sem <- error_95 <- NULL
  error_bars <- match.arg(error_bars)
  if (.is_na(colors)) {
    # colors <- COLORS7
    colors <- c("#2EBAED", "#000000", "#DE1C14", "#E98C7B", "#D4D2D2", "#ADCC54", 
                "#F0D0CE")
  }
  if (length(colors) != 7) {
    stop("Colors parameter: supply color vector with length 7")
  }
  if (CT == FALSE) {
    type_occurrences <- type_occurrences[, seq_len(6)]
  }
  else {
    type_occurrences <- type_occurrences[, c(1, 2, 8, 7, 
                                             4, 5, 6)]
  }
  if (.is_na(by)) {
    by <- "all"
  }
  tb_per_sample <- type_occurrences %>% tibble::rownames_to_column("sample") %>% 
    dplyr::mutate(by = by) %>% tidyr::pivot_longer(c(-sample, 
                                                     -by), names_to = "variable", values_to = "nmuts") %>% 
    dplyr::group_by(sample) %>% dplyr::mutate(value = nmuts) %>% 
    dplyr::ungroup() %>% dplyr::mutate(sub_type = stringr::str_remove(variable, 
                                                                      " .*"), variable = factor(variable, levels = unique(variable)))
  tb <- tb_per_sample %>% dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, variable) %>% dplyr::summarise(sub_type = sub_type[[1]], 
                                                       mean = mean(value), stdev = stats::sd(value), total_individuals = sum(value), group_number = length(value), type_sum = mean * group_number,
                                                       total_mutations = sum(nmuts)) %>% dplyr::mutate(total_individuals = sum(total_individuals), 
                                                                                                       total_mutations = sum(total_mutations), total_mutations_per_cell = round(total_individuals / group_number, digits = 0)) %>% dplyr::mutate(sem = stdev/sqrt(total_individuals), 
                                                                                                                                                                                                                                                 error_95 = ifelse(total_individuals > 1, qt(0.975, df = total_individuals - 
                                                                                                                                                                                                                                                                                               1) * sem, NA)) %>% dplyr::ungroup() %>% dplyr::mutate(total_mutations = prettyNum(total_mutations, 
                                                                                                                                                                                                                                                                                                                                                                                 big.mark = ","), total_mutations = paste("Net change (IHD - Control) = ", 
                                                                                                                                                                                                                                                                                                                                                                                                                          total_mutations_per_cell), error_pos = mean)
  # view(tb)
  if (CT == FALSE) {
    colors <- colors[c(1, 2, 3, 5, 6, 7)]
  }
  else {
    CpG <- which(tb$variable == "C>T at CpG")
    other <- which(tb$variable == "C>T other")
    tb$error_pos[CpG] <- tb$error_pos[other] + tb$error_pos[CpG]
    CpG <- which(tb_per_sample$variable == "C>T at CpG")
    other <- which(tb_per_sample$variable == "C>T other")
    tb_per_sample$value[CpG] <- tb_per_sample$value[other] + 
      tb_per_sample$value[CpG]
  }
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.9
    spacing <- 0.5
  }
  plot <- ggplot(data = tb, aes(x = sub_type, y = mean, fill = variable, 
                                group = sub_type, width = width)) + geom_bar(stat = "identity") + 
    # geom_text(aes(label=round(mean, digits = 0)), vjust=-0.25) +
    scale_fill_manual(values = colors, name = "Point mutation type") + 
    theme_bw() + xlab("") + ylab("number of sSNV") + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), 
          panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
                                                                       "lines"))
  if (indv_points == TRUE) {
    tb_per_sample <- dplyr::left_join(tb_per_sample, tb[, 
                                                        c("by", "variable", "total_mutations")], by = c("by", 
                                                                                                        "variable"))
    plot <- plot + geom_jitter(data = tb_per_sample, aes(y = value), 
                               height = 0, width = 0.3, shape = 21, colour = "grey23")
  }
  if (sum(is.na(tb$stdev)) > 0 & error_bars != "none") {
    warning("No error bars can be plotted, because there is only one sample per mutation spectrum.\n              Use the argument: `error_bars = 'none'`, if you want to avoid this warning.", 
            call. = FALSE)
  }
  else {
    if (error_bars == "stdev") {
      plot <- plot + geom_errorbar(aes(ymin = error_pos - 
                                         stdev, ymax = error_pos + stdev), width = 0.2)
    }
    else if (error_bars == "95%_CI") {
      plot <- plot + geom_errorbar(aes(ymin = error_pos - 
                                         error_95, ymax = error_pos + error_95), width = 0.2)
    }
    else if (error_bars == "SEM") {
      plot <- plot + geom_errorbar(aes(ymin = error_pos - 
                                         sem, ymax = error_pos + sem), width = 0.2)
    }
  }
  if (length(by) == 1) {
    plot <- plot + facet_wrap(~total_mutations)
  }
  else {
    plot <- plot + facet_wrap(factor(by, c("Control", "Disease")) ~ total_mutations)
    # plot <- plot + facet_wrap(by ~ total_mutations)
  }
  if (legend == FALSE) {
    plot <- plot + guides(fill = "none")
  }
  return(plot)
}

.is_na <- function (x) 
{
  purrr::is_scalar_vector(x) && is.na(x)
}

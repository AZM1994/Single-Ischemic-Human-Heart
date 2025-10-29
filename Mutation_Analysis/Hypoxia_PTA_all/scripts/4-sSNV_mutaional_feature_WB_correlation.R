library(tibble)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)

## WB data
WB_data <- tribble(
  ~donor, ~MSH2, ~MSH6, ~MLH3, ~MLH1, ~Ku80, ~XRCC1, ~OGG1, ~HIF1A, ~C_A_Substitution, ~Total_sSNV_burden,
  1863, 5153, 6929, 9071, 7902, 13085, 8140, 3746, 4233, 22.0, 123.0,
  1039, 3923, 6812, 9896, 8277, 15419, 7177, 8902, 9005, 67.8, 350.9,
  1940, 2769, 3993, 8782, 8078, 14807, 5865, 10660, 10140, 27.3, 226.7,
  5919,  521, 1736, 6400, 6159, 14147, 6099, NA,    NA,   74.3, 405.5,
  1363,  127, 1258, 3275, 3745, 11804, 3785, 16904, 17281, 83.1, 440.7,
  1673,   99, 1390, 2659, 2440, 7005, 1380, 19960, 19432, 213.6, 639.3,
  604,  126, 1317, 2393, 2135, 4871,  384, 19354, 16477, 172.9, 775.1,
  1743,  125, 1230, 2370, 2425, 3446,  363, NA,    NA,  438.7, 3012.9,
  1113, 1123, 1096, 2244, 2429, 2566,  243, 21521, 17750, 231.7, 808.6
)

## mutation features
mutation_features <- tribble(
  ~donor, ~SBS1, ~SBS4, ~SBS5, ~SBS8, ~SBS18, ~SBS30, ~SBS32, ~SBS36, ~SBS40a, ~SBS44, ~SBS89, ~SBS92,
  1863,  5.3,  0.2,  52.0,  3.4,  2.3,  0.3,  3.5,  0.0,  0.4,  0.0,  5.4,  0.0,
  1039, 38.0,  0.0, 246.1, 38.1,  0.0,  0.6,  0.6,  0.0,  0.0,  0.0,  0.0,  0.0,
  1940,  5.6,  0.8, 133.3,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  5919, 38.8,  0.0, 228.3,  6.6,  0.0,  2.6,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  1363, 29.5,  0.0, 114.2, 25.5,  0.0, 39.7, 45.4,  3.5, 17.5, 15.3, 19.2,  1.3,
  1673, 26.2,  51.0, 443.2, 58.2,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
  604,  10.8,  59.1, 188.3, 72.7, 13.3, 46.9, 18.9, 58.6,122.0, 38.2, 17.0,  0.0,
  1743, 175.6, 124.9,1622.1, 26.8,  1.4,131.7, 67.2,  0.0,162.7, 63.2, 41.0,138.7,
  1113, 27.1,  60.8, 590.2, 39.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0
)

WB_mutation_features <- WB_data %>% left_join(mutation_features, by = "donor")

# Define predictors and outcomes
wb_vars <- c("MSH2","MSH6","MLH3","MLH1","Ku80","XRCC1","OGG1","HIF1A")
# mutation_features_vars <- c("SBS1","SBS4","SBS5","SBS8","SBS18","SBS30","SBS32","SBS36","SBS40a","SBS44","SBS89","SBS92","C_A_Substitution","Total_sSNV_burden")
# mutation_features_vars <- c("SBS4","SBS8","SBS18","SBS30","SBS32","SBS36","SBS40a","SBS44","C_A_Substitution","Total_sSNV_burden")
mutation_features_vars <- c("SBS4","SBS8","SBS30","SBS44","C_A_Substitution","Total_sSNV_burden")

## Function for Spearman correlation
get_corr <- function(x, y, data) {
  df <- data %>% dplyr::select(all_of(c(x, y))) %>% na.omit()
  if(nrow(df) > 2) {
    res <- cor.test(df[[x]], df[[y]], method = "spearman")
    tibble(
      predictor = x,
      outcome = y,
      rho = unname(res$estimate),
      pval = res$p.value,
      n = nrow(df)
    )
  } else {
    tibble(predictor = x, outcome = y, rho = NA, pval = NA, n = nrow(df))
  }
}

## Run all correlations
results <- map_dfr(wb_vars, function(wb) {
  map_dfr(mutation_features_vars, function(out) {
    get_corr(wb, out, WB_mutation_features)
  })
}) %>% 
  mutate(padj = p.adjust(pval, method = "fdr"))
write.csv(results, paste0(table_dir, "/4-WB_mutational_features_correlation_test.csv"), row.names = F)

## Define groups
row_order <- c("OGG1","HIF1A", "MSH2","MSH6","MLH1","MLH3", "XRCC1","Ku80")
row_group <- c("Hypoxia/oxidative stress","Hypoxia/oxidative stress","MMR","MMR","MMR","MMR","BER","BER")
row_groups_df <- tibble(predictor = row_order, group = row_group)

col_order <- c("C_A_Substitution","Total_sSNV_burden","SBS4","SBS8","SBS30","SBS44")

plot_df <- results %>% 
  mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "*", ""), 
         predictor = factor(predictor, levels = row_order), 
         outcome = factor(outcome, levels = col_order)) %>% 
  left_join(row_groups_df, by = "predictor") %>% 
  mutate(group = factor(group, levels = c("Hypoxia/oxidative stress","MMR","BER")))

## Heatmap
p_WB_mutational_features <- ggplot(plot_df, aes(x = outcome, y = predictor, fill = rho)) + 
  geom_tile(color = "grey80") + geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5, size = 8) + 
  scale_fill_gradient2(low = "#246aae", mid = "white", high = "#870b24", midpoint = 0, limits = c(-1,1), name = expression("Spearman " * rho)) + 
  facet_grid(group ~ ., scales = "free", space = "free") + theme_minimal(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold")) + 
  labs(x = "Mutational features", y = "DNA repair and hypoxia proteins")
ggsave(paste0(main_figure_dir, "/4-WB_mutational_features.pdf"), plot = p_WB_mutational_features, width = 9, height = 6, dpi = 600)

#############################################################################
################ compare the original and revised SCAN2 call ################
#############################################################################

##### read revised snv and indel
revised_snv_file_list <- list.files(path = paste0(project_dir, "/new_mut_burden_table/revised/snv"), pattern = ".csv", full.names = TRUE)
revised_snv_summary_df <- rbindlist(sapply(revised_snv_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
revised_snv_summary_df <- revised_snv_summary_df[revised_snv_summary_df$V1 == 2]

revised_snv_summary_df$cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", revised_snv_summary_df$cell_ID))
revised_snv_summary_df$cell_ID <- gsub('snv_gatk_table_', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('_revised', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('1_', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('-2n', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('M-', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('_2n', "", revised_snv_summary_df$cell_ID)
revised_snv_summary_df$cell_ID <- gsub('CM_', "", revised_snv_summary_df$cell_ID)

revised_indel_file_list <- list.files(path = paste0(project_dir, "/new_mut_burden_table/revised/indel"), pattern = ".csv", full.names = TRUE)
revised_indel_summary_df <- rbindlist(sapply(revised_indel_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
revised_indel_summary_df <- revised_indel_summary_df[revised_indel_summary_df$V1 == 2]

revised_indel_summary_df$cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", revised_indel_summary_df$cell_ID))
revised_indel_summary_df$cell_ID <- gsub('indel_gatk_table_', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('_revised', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('1_', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('-2n', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('M-', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('_2n', "", revised_indel_summary_df$cell_ID)
revised_indel_summary_df$cell_ID <- gsub('CM_', "", revised_indel_summary_df$cell_ID)

##### read original snv and indel
scan2_snv_file_list <- list.files(path = paste0(project_dir, "/new_mut_burden_table/scan2/snv"), pattern = ".csv", full.names = TRUE)
scan2_snv_summary_df <- rbindlist(sapply(scan2_snv_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
scan2_snv_summary_df <- scan2_snv_summary_df[scan2_snv_summary_df$V1 == 2]

scan2_snv_summary_df$cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", scan2_snv_summary_df$cell_ID))
scan2_snv_summary_df$cell_ID <- gsub('snv_gatk_table_', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('_scan2', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('1_', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('-2n', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('M-', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('_2n', "", scan2_snv_summary_df$cell_ID)
scan2_snv_summary_df$cell_ID <- gsub('CM_', "", scan2_snv_summary_df$cell_ID)

scan2_indel_file_list <- list.files(path = paste0(project_dir, "/new_mut_burden_table/scan2/indel"), pattern = ".csv", full.names = TRUE)
scan2_indel_summary_df <- rbindlist(sapply(scan2_indel_file_list, fread, simplify = FALSE), idcol = 'cell_ID')
scan2_indel_summary_df <- scan2_indel_summary_df[scan2_indel_summary_df$V1 == 2]

scan2_indel_summary_df$cell_ID <- sub(".*\\\\([^.]*).*", "\\1", gsub(".*/(.+).csv*", "\\1", scan2_indel_summary_df$cell_ID))
scan2_indel_summary_df$cell_ID <- gsub('indel_gatk_table_', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('_scan2', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('1_', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('-2n', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('M-', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('_2n', "", scan2_indel_summary_df$cell_ID)
scan2_indel_summary_df$cell_ID <- gsub('CM_', "", scan2_indel_summary_df$cell_ID)

##### put snv and indel together
redo_scan2_summary <- cbind(scan2_snv_summary_df, scan2_indel_summary_df[,-(1:2)])
colnames(redo_scan2_summary) <- colnames(SCAN2_df)
# compare_01 <- redo_scan2_summary[,3:10] - SCAN2_df[,3:10]
# compare_02 <- redo_scan2_summary[,12:19] - SCAN2_df[,12:19]

revised_summary <- cbind(revised_snv_summary_df, revised_indel_summary_df[,-(1:2)])
colnames(revised_summary) <- colnames(SCAN2_df)
# compare_03 <- revised_summary[,3:10] - SCAN2_df[,3:10]
# compare_04 <- revised_summary[,12:19] - SCAN2_df[,12:19]

redo_scan2_summary$type <- 'redo_scan2'
revised_summary$type <- 'revised'
SCAN2_df$type <- 'ori_scan2'

# all_together <- rbind(SCAN2_df, redo_scan2_summary, revised_summary)
all_together <- rbind(revised_summary, redo_scan2_summary)

pdf(paste0(project_dir, "/new_mut_burden_table", "/compare_revised_old_mut_burden.pdf"), width = 20, height = 10)
p1 <- ggplot(all_together, aes(x = factor(cell_ID, level = sample_name_list), y = snv.rate.per.gb, color = type, shape = type)) +
  geom_point(size = 3, na.rm = FALSE) +
  scale_y_break(c(1500, 3300)) + ylim(0, 3500) +
  labs(x = "Sample ID", y = "sSNV rate per GB") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))
p1

p2 <- ggplot(all_together, aes(x = factor(cell_ID, level = sample_name_list), y = indel.rate.per.gb, color = type, shape = type)) +
  geom_point(size = 3, na.rm = FALSE) +
  labs(x = "Sample ID", y = "sindel rate per GB") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 0.0, hjust = 0.0))
p2
dev.off()
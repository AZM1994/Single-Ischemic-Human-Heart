# Rscript /n/no_backup2/bch/lee/zheming/heart_hypoxia/scan2_analysis/Scan2_R_Analysis/load_scan2_rda.R B4 IPSC_Hypoxia_PTA
args <- commandArgs(TRUE)
verbose  = TRUE
cell_id  = args[1]
case_dir = args[2]

load(sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/%s/call_mutations/%s/scan2_object.rda",
             case_dir,cell_id))

sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/%s/call_mutations/%s/scan2_object.rda",
        case_dir,cell_id)
sink(file = sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/scan2_output/scan2_output_%s.txt",cell_id))
results
sink(file = NULL)

write.csv(results@mutburden, sprintf("/n/no_backup2/bch/lee/zheming/heart_hypoxia/IPSC_Cases_scan2_analysis/Scan2_R_Analysis/all_burden_table/%s.csv",cell_id), row.names=TRUE)
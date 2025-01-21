# Data Folder Structure

- `data/`: Contains all the data needed for SigNet
  - `data.xlsx`: [COSMIC v3.1 Signatures](https://cancer.sanger.ac.uk/signatures/downloads/)
  - `data_v2.xlsx`: [COSMIC v2.0 Signatures](https://cancer.sanger.ac.uk/signatures/downloads/)
  - `mutation_type_order.xsls`: Mapping to correctly sort mutation types
  - `datasets/` 
    - `real_input.csv`, `real_label.csv`: Input and labels used to benchmark all methods
    - `detector/`: Contains all datasets needed to train the Detector module
    - `finetuner/`: Contains all datasets needed to train the Finetuner module
    - `errorfinder/`: Contains all datasets needed to train the Errorfinder module
  - `real_data/`:
    - `3mer_WG_hg37.txt`: Whole genome human abundances with reference genome GRCh37
    - `3mer_WG_hg38.txt`: Whole genome human abundances with reference genome GRCh38
    - `abundancies_trinucleotides.txt`: Whole exon human abundances
    - `PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv`: [Real total counts](https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples)
    - `sigprofiler_not_norm_PCAWG.csv`: Real weights used to create the training sets (normalizing the real total counts)
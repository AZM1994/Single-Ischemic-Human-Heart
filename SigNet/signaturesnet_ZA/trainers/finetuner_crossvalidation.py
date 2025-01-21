import os
import sys

import numpy as np
import pandas as pd
import torch

from signaturesnet.utilities.oversampler import CancerTypeOverSampler
from signaturesnet.utilities.metrics import get_classification_metrics
from signaturesnet.utilities.io import csv_to_pandas, read_config, read_signatures
from signaturesnet.utilities.data_partitions import DataPartitions
from signaturesnet.utilities.data_generator import DataGenerator
from signaturesnet.trainers import train_finetuner
from signaturesnet.modules import CombinedFinetuner
from signaturesnet.models import Baseline

def partition_dataset(df, n_partitions):
    df["index"] = list(range(len(df)))
    df["fold"] = df["index"]%n_partitions
    return [df[df["fold"] == k].drop(["index", "fold"], axis=1) for k in range(n_partitions)]

def read_data_and_partitions(k):
    '''
    Read real data and generate k partitions to apply k-fold crossvalidation.
    '''
    dev = "cuda" if torch.cuda.is_available() else "cpu"

    real_data = csv_to_pandas("../../data/real_data/sigprofiler_not_norm_PCAWG.csv",
                              device=dev, header=0, index_col=0,
                              type_df="../../data/real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv")

    # Shuffle samples inside the same cancer type
    real_data = real_data.groupby('cancer_type').sample(frac=1, random_state=0)

    # Partition the data in k groups preserving proportion of samples in each cancer type
    lst = partition_dataset(real_data, n_partitions=k)

    # Remove cancer column, normalize each row and add missing signatuers with 0s at the end
    lst_weights = [torch.tensor(df.values[:, :-1], dtype=torch.float) for df in lst]
    lst_weights = [df/torch.sum(df, axis=1).reshape(-1, 1) for df in lst_weights]
    lst_weights = [torch.cat([df, torch.zeros(df.size(0), 7).to(df)], dim=1) for df in lst_weights]

    # Cancer type column
    lst_ctype = [torch.tensor(df.values[:, -1], dtype=torch.float) for df in lst]
    return lst_weights, lst_ctype


if __name__ == "__main__":

    k = 10
    assert len(sys.argv) == 2
    network_type = str(sys.argv[1])
    assert network_type in ['low', 'large']
    finetuner_config_path = "../configs/finetuner/finetuner_" + network_type + ".yaml"

    # Create partitions
    lst_weights, lst_ctype = read_data_and_partitions(k)

    # Create inputs associated to the labels:
    signatures = read_signatures(
        "../../data/data.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
    data_generator = DataGenerator(signatures=signatures,
                                   seed=0,
                                   shuffle=True)
    
    N_oversample_list = [5, 10, 20, 30]
    
    for N_oversample in N_oversample_list:
        # Oversample each set to have the same number of samples for each cancer type
        oversampled_weights = []
        for i in range(k):
            oversampler = CancerTypeOverSampler(lst_weights[i], lst_ctype[i])
            oversampled_weights.append(oversampler.get_N_oversampled_set(N_oversample))
        
        # Loop through the partitions
        for i in range(k):
            current_test = i
            current_val = (i+1)%k
            current_train = [j for j in range(k) if j != current_test and j != current_val]

            # Create train, val and test weight sets
            test_weights = oversampled_weights[current_test]
            val_weights = oversampled_weights[current_val]
            train_weights = torch.cat([oversampled_weights[j] for j in current_train], axis=0)

            # Create pairs input-label
            print("Creating train, val and test data")
            train_input, train_label = data_generator.make_input(train_weights, "train", network_type, normalize=True)
            val_input, val_label = data_generator.make_input(val_weights, "val", network_type, normalize=True)
            test_input, test_label = data_generator.make_input(test_weights, "test", network_type, normalize=True)

            # Run Baseline
            print("Running Baseline")
            sf = Baseline(signatures)
            train_baseline = sf.get_weights_batch(train_input, n_workers=10)
            val_baseline = sf.get_weights_batch(val_input, n_workers=10)
            test_baseline = sf.get_weights_batch(test_input, n_workers=10)

            # Create DataPartitions
            train_data = DataPartitions(inputs=train_input,
                                        prev_guess=train_baseline,
                                        labels=train_label)
            val_data = DataPartitions(inputs=val_input,
                                    prev_guess=val_baseline,
                                    labels=val_label)

            # Read config
            config = read_config(path=finetuner_config_path)
            config["wandb_project_id"] = "crossval_finetuner"
            config["models_dir"] = "../../trained_models/crossval_oversample"
            config["model_id"] = "finetuner_" + network_type + "_crossval_" + str(i) + "_" + str(N_oversample)
            train_data.to(config["device"])
            val_data.to(config["device"])

            # Train finetuner
            print("Training Finetuner")
            score = train_finetuner(config=config, train_data=train_data, val_data=val_data)
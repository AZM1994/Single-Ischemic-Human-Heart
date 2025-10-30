import logging
import os
import sys

import torch
import matplotlib.pyplot as plt
import numpy as np

from signaturesnet.utilities.data_generator import DataGenerator
from signaturesnet.utilities.io import read_model, read_data_generator, read_signatures, tensor_to_csv

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Usage: python generate_data_realistic v3"
    cosmic_version = str(sys.argv[1])

    if cosmic_version == 'v3':
        generator_model_path = "../../trained_models/exp_not_norm/generator"
        experiment_id = "exp_not_norm"
        signatures = read_signatures("../../data/data.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
    elif cosmic_version == 'v2':
        generator_model_path = "../../trained_models/exp_generator_v2/generator"
        experiment_id = "exp_generator_v2"
        signatures = read_signatures("../../data/data_v2.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
    else:
        raise NotImplementedError

    num_samples = 5000

    data_generator = DataGenerator(signatures=signatures,
                                   seed=0,
                                   shuffle=True)

    datasets = [("low", "train"), ("low", "val"), ("large", "train"), ("large", "val"), ("large", "test")]
    for large_low, split in datasets:
        logging.info("Generating dataset %s, %s"%(large_low, split))
        n_samples = int(num_samples/100) if split in ["val", "test"] else num_samples
        train_input, train_labels = data_generator.make_realistic_set(generator_model_path=generator_model_path,
                                                                      large_low=large_low,
                                                                      split=split,
                                                                      num_samples=n_samples,
                                                                      std=1.5)
        logging.info("Dataset created")
        tensor_to_csv(train_input, "../../data/%s/sd1.5/%s_generator_%s_input.csv"%(experiment_id, split, large_low))
        logging.info("Saved inputs")
        tensor_to_csv(train_labels, "../../data/%s/sd1.5/%s_generator_%s_label.csv"%(experiment_id, split, large_low))
        logging.info("Saved labels")
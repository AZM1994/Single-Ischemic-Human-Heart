import os
import sys

import numpy as np
import torch

from signaturesnet.utilities.io import read_signatures, tensor_to_csv, csv_to_tensor
from signaturesnet.utilities.data_generator import DataGenerator
from signaturesnet.models import Baseline

def shuffle(inputs, labels, num_mut):
    indexes = torch.randperm(inputs.shape[0])
    return inputs[indexes, ...], labels[indexes, ...], num_mut[indexes, ...]

if __name__ == "__main__":

    assert len(sys.argv) == 2, "Usage: python generate_classifier_data.py v3"
    cosmic_version = str(sys.argv[1])

    if cosmic_version == 'v3':
        experiment_id = "exp_errorfinder_final"
        signatures = read_signatures("../../data/data.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
    elif cosmic_version == 'v2':
        experiment_id = "exp_errorfinder_final"
        signatures = read_signatures("../../data/data_v2.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
    else:
        raise NotImplementedError

    data_folder = "../../data"

    data_generator = DataGenerator(signatures=signatures,
                                    seed=0,
                                    shuffle=True)

    real_data = csv_to_tensor(data_folder + "/real_data/sigprofiler_not_norm_PCAWG.csv", header=0, index_col=0)
    real_weights = torch.cat([real_data, torch.zeros(real_data.size(0), 7).to(real_data)], dim=1)

    for s, ll in [("train", "low"), ("train", "large"), ("val", "low"), ("val", "large"), ("test", None)]:
        inputs, labels = data_generator.make_input(labels=real_weights, split=s, large_low=ll)
        baseline = Baseline(signatures=signatures).get_weights_batch(inputs)
        if ll is not None:
            tensor_to_csv(inputs, data_folder + '/' + experiment_id + "/%s_%s_input.csv"%(s, ll))
            tensor_to_csv(labels, data_folder + '/' + experiment_id + "/%s_%s_label.csv"%(s, ll))
            tensor_to_csv(baseline, data_folder + '/' + experiment_id + "/%s_%s_baseline.csv"%(s, ll))
        else:
            tensor_to_csv(inputs, data_folder + '/' + experiment_id + "/%s_input.csv"%(s))
            tensor_to_csv(labels, data_folder + '/' + experiment_id + "/%s_label.csv"%(s))
            tensor_to_csv(baseline, data_folder + '/' + experiment_id + "/%s_baseline.csv"%(s))

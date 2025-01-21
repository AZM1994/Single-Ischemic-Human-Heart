import os
import sys

import numpy as np
import torch

from signaturesnet.utilities.io import read_signatures, tensor_to_csv, csv_to_tensor
from signaturesnet.utilities.data_generator import DataGenerator
from signaturesnet import DATA

def shuffle(inputs, labels, num_mut):
    indexes = torch.randperm(inputs.shape[0])
    return inputs[indexes, ...], labels[indexes, ...], num_mut[indexes, ...]

if __name__ == "__main__":

    assert len(sys.argv) == 2, "Usage: python generate_classifier_data.py v3"
    cosmic_version = str(sys.argv[1])

    data_folder = DATA

    if cosmic_version == 'v3':
        experiment_id = "exp_classifier_superlow"
        signatures = read_signatures(data_folder + "/data.xlsx", mutation_type_order=data_folder + "/mutation_type_order.xlsx")
    elif cosmic_version == 'v2':
        experiment_id = "exp_classifier_superlow"
        signatures = read_signatures(data_folder + "/data_v2.xlsx", mutation_type_order=data_folder + "/mutation_type_order.xlsx")
    else:
        raise NotImplementedError

    data_generator = DataGenerator(signatures=signatures,
                                   seed=0,
                                   shuffle=True)

    # Read
    real_data = csv_to_tensor(data_folder + "/real_data/sigprofiler_not_norm_PCAWG.csv", header=0, index_col=0)
    real_weights = torch.cat([real_data, torch.zeros(real_data.size(0), 7).to(real_data)], dim=1)
    real_weights = real_weights.repeat(2, 1)
    realistic_input, realistic_weights = data_generator.make_input(labels=real_weights,
                                                                   split="train",
                                                                   large_low="superlow")
    random_input, random_weights = data_generator.make_random_set(split='train',
                                                    large_low='superlow',
                                                    num_samples=int(real_weights.shape[0]),
                                                    min_n_signatures=1,
                                                    max_n_signatures=10)
    realistic_nummut = realistic_weights[:, -1].view(-1, 1)
    random_nummut = random_weights[:, -1].view(-1, 1)

    # Label realistic data as 1, random as 0
    realistic_labels = torch.ones((realistic_input.shape[0], 1)).to(torch.float).view(-1, 1)
    random_labels = torch.zeros((random_input.shape[0], 1)).to(torch.float).view(-1, 1)

    # Concatenate and shuffle datasets
    inputs = torch.cat((realistic_input, random_input))
    labels = torch.cat((realistic_labels, random_labels))
    num_mut = torch.cat((realistic_nummut, random_nummut))
    inputs, labels, num_mut = shuffle(inputs, labels, num_mut)
    
    # Partition datasets
    n = inputs.shape[0]
    train_classification_inputs = inputs[:int(n*0.7), ...]
    train_classification_labels = labels[:int(n*0.7), ...]
    train_classification_numut = num_mut[:int(n*0.7), ...]

    val_classification_inputs = inputs[int(n*0.7):int(n*0.85), ...]
    val_classification_labels = labels[int(n*0.7):int(n*0.85), ...]
    val_classification_numut = num_mut[int(n*0.7):int(n*0.85), ...]

    test_classification_inputs = inputs[int(n*0.85):, ...]
    test_classification_labels = labels[int(n*0.85):, ...]
    test_classification_numut = num_mut[int(n*0.85):, ...]

    # store everything
    tensor_to_csv(train_classification_inputs, data_folder + '/' + experiment_id + "/classifier/train_input.csv")
    tensor_to_csv(train_classification_labels, data_folder + '/' + experiment_id + "/classifier/train_label.csv")
    tensor_to_csv(train_classification_numut, data_folder + '/' + experiment_id + "/classifier/train_num_mut.csv")

    tensor_to_csv(val_classification_inputs, data_folder + '/' + experiment_id + "/classifier/val_input.csv")
    tensor_to_csv(val_classification_labels, data_folder + '/' + experiment_id + "/classifier/val_label.csv")
    tensor_to_csv(val_classification_numut, data_folder + '/' + experiment_id + "/classifier/val_num_mut.csv")

    tensor_to_csv(test_classification_inputs, data_folder + '/' + experiment_id + "/classifier/test_input.csv")
    tensor_to_csv(test_classification_labels, data_folder + '/' + experiment_id + "/classifier/test_label.csv")
    tensor_to_csv(test_classification_numut, data_folder + '/' + experiment_id + "/classifier/test_num_mut.csv")



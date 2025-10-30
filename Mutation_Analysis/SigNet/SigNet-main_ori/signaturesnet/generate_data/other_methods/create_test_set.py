import os
import sys

import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utilities.io import csv_to_tensor, read_signatures, tensor_to_csv
from utilities.data_generator import DataGenerator

experiment_id = "exp_all"

# Creat test set:
real_weights = csv_to_tensor("../../../data/real_data/sigprofiler_not_norm_PCAWG.csv", device='cpu', header=0, index_col=0)
real_weights = torch.cat([real_weights, torch.zeros(real_weights.size(0), 7).to(real_weights)], dim=1) 

signatures = read_signatures(
    "../../../data/data.xlsx", mutation_type_order="../../../data/mutation_type_order.xlsx")
datagenerator = DataGenerator(signatures)
test_input, test_label = datagenerator.make_input(labels = real_weights, set = "test", large_low = "large", normalize=True)

tensor_to_csv(test_input, "../../../data/%s/test_input.csv"%experiment_id)
tensor_to_csv(test_label, "../../../data/%s/test_label.csv"%experiment_id)
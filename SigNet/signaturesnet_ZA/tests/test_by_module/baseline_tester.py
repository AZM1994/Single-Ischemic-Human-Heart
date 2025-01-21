import os
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from models import Baseline
from utilities.io import csv_to_tensor, read_signatures
from utilities.plotting import plot_all_metrics_vs_mutations,\
                               plot_metric_vs_mutations,\
                               plot_interval_metrics_vs_mutations,\
                               plot_interval_performance
from modules import CombinedFinetuner, CombinedErrorfinder

# Load data
data_path = "../../../data/exp_good/test_realistic/"

inputs = csv_to_tensor(data_path + "test_realistic_input.csv", device='cpu')
# baselines = csv_to_tensor(data_path + "test_realistic_baseline.csv", device='cpu')
baselines = csv_to_tensor(data_path + "test_realistic_label.csv", device='cpu')[:-1, :-1]
labels = csv_to_tensor(data_path + "test_realistic_label.csv", device='cpu')[:-1, ...]
num_mut = labels[:, -1].unsqueeze(1)

# from scipy import stats
# import torch
# print(stats.describe(torch.sum(labels[..., :-1], axis=1)))
# print(labels[..., :-1]))

# import matplotlib.pyplot as plt
# plt.hist(torch.sum(baselines, axis=1), density=True, bins=30)
# plt.show()

print(baselines.shape)
print(labels.shape)

signatures = read_signatures(file="../../../data/data.xlsx",
                             mutation_type_order="../../../data/mutation_type_order.xlsx")

# Plot results
list_of_methods = ['baseline']
list_of_guesses = [baselines]
# plot_all_metrics_vs_mutations( list_of_methods, list_of_guesses, labels, '')
plot_metric_vs_mutations(list_of_metrics=["accuracy %", "reconstruction_error"],
                         list_of_methods=list_of_methods,
                         list_of_guesses=list_of_guesses,
                         label=labels,
                         show=True,
                         signatures=signatures,
                         mutation_distributions=inputs)

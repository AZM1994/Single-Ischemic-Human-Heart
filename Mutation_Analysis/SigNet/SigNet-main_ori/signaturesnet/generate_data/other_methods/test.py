import matplotlib.pyplot as plt
import pandas as pd
import torch
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.signaturesnet import SigNet
from utilities.io import csv_to_tensor, read_signatures
from utilities.plotting import plot_reconstruction


input_file_path = "../../data/real_data/PCAWG_data.csv"
opportunities = "genome"
output_path = "../../tmp/real_data_final_3" 
plot_figs = False

signatures = read_signatures("../../data/data.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")


# NOT NORMALIZED DATA:

labels = csv_to_tensor("../../data/real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv",
                              device="cpu", header=0, index_col=0)
labels = torch.cat([labels, torch.zeros(labels.size(0), 7).to(labels)], dim=1)

input_file = pd.read_csv(input_file_path, header=0, index_col=0)
mutation_data = torch.tensor(input_file.values, dtype=torch.float)

# mutation_data/torch.sum(mutation_data, axis=1).reshape(-1, 1).detach().numpy()
plot_reconstruction(mutation_data.detach().numpy(), labels, signatures, [100,101,200], '')


# NORMALIZED DATA:
labels = csv_to_tensor("../../data/real_data/sigprofiler_normalized_PCAWG.csv",
                              device="cpu", header=0, index_col=0)
labels = labels/torch.sum(labels, axis=1).reshape(-1, 1)
labels = torch.cat([labels, torch.zeros(labels.size(0), 7).to(labels)], dim=1)

signet = SigNet(opportunities_name_or_path=opportunities, signatures_path="../../data/data.xlsx")

input_file = pd.read_csv(input_file_path, header=0, index_col=0)
mutation_data = torch.tensor(input_file.values, dtype=torch.float)
weight_guess, upper_bound, lower_bound, classification, normalized_input = signet(mutation_vec=mutation_data)

plot_reconstruction(normalized_input, labels, signatures, [100,101,200], '')



# EXACT DATA: IT WORKS.
# labels = csv_to_tensor("../../data/real_data/exact_real_data_labels.csv")
# labels = labels[:,:-1]
# input_file = pd.read_csv("../../data/real_data/exact_real_data_input.csv", header=None, index_col=None)
# mutation_data = torch.tensor(input_file.values, dtype=torch.float)
# plot_reconstruction(mutation_data.detach().numpy(), labels, signatures, [0,1,2], '')

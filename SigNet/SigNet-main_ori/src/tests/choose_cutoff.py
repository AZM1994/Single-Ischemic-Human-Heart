from email.mime import base
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from models.baseline import Baseline
from modules.combined_finetuner import CombinedFinetuner
from utilities.io import csv_to_tensor, read_signatures
from utilities.metrics import get_MSE


 # This is for the test set
# input_batch = csv_to_tensor("../../data/exp_all/test_input.csv")
# label_batch = csv_to_tensor("../../data/exp_all/test_label.csv")

 # This is for the training set
input_batch = csv_to_tensor("../../data/exp_all/errorfinder/train_low_input.csv")
label_batch = csv_to_tensor("../../data/exp_all/errorfinder/train_low_label.csv")

signatures = read_signatures("../../data/data.xlsx")
# Run Baseline
print("Running Baseline")
sf = Baseline(signatures)
test_baseline = sf.get_weights_batch(input_batch, n_workers=4)

models_path = "../../trained_models/exp_final"
finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "/finetuner_low",
                                large_mum_mut_dir=models_path + "/finetuner_large")
finetuner_guess = finetuner(mutation_dist=input_batch,
                            baseline_guess=test_baseline,
                            num_mut=label_batch[:, -1].view(-1, 1))

cutoffs = np.linspace(0,0.1,10)
cutoffs[0] = 0.001  # This is the minimum cutoff that we put in the finetuner.
num_muts = [15,50,100,250,500,1000,5e3,1e4]                # This is for the training set
# num_muts = np.unique(label_batch[:,-1].detach().numpy())              # This is for the test set
mse_cutoffs = np.zeros((len(cutoffs), len(num_muts)))
fp_cutoffs = np.zeros((len(cutoffs), len(num_muts)))

for j, cutoff in enumerate(cutoffs):
    for i, num_mut in enumerate(num_muts):
        try:
            indexes = label_batch[:, -1] < num_muts[i+1] and label_batch[:, -1] > num_mut       # This is for the training set
        except:
            indexes = label_batch[:, -1] > num_mut                                              # This is for the training set
        # indexes = label_batch[:, -1] == num_mut                               # This is for the test set

        guess_cutoff = finetuner_guess.clone().detach()
        guess_cutoff[guess_cutoff<cutoff] = 0
        mse_cutoffs[j, i] = get_MSE(guess_cutoff[indexes, :], label_batch[indexes, :-1])
        fp_cutoffs[j, i] = torch.sum(guess_cutoff[indexes, :][label_batch[indexes, :-1]<1e-4]>0)

plt.plot(np.log10(num_muts),np.transpose(mse_cutoffs), marker='o',linewidth=0.5)
plt.xlabel('log(N)')
plt.ylabel('MSE')
plt.legend(cutoffs)
plt.show()

plt.plot(np.log10(num_muts),np.transpose(fp_cutoffs), marker='o',linewidth=0.5)
plt.xlabel('log(N)')
plt.ylabel('FP')
plt.legend(cutoffs)
plt.show()
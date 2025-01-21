
from cmath import nanj
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sn
import torch 

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from models import Baseline
from modules import CombinedFinetuner
from tests.real_data import read_finetuner
from trainers.finetuner_crossvalidation import read_data_and_partitions
from utilities.data_generator import DataGenerator
from utilities.io import csv_to_pandas, csv_to_tensor, read_signatures, sort_signatures
from utilities.metrics import get_MSE
from utilities.oversampler import CancerTypeOverSampler
from utilities.plotting import get_correlation_matrix, plot_correlation_matrix


def bootstrap_std(set, bootstrap_size, N_bootstrap, signatures):
    corr_list = []
    for _ in range(N_bootstrap):
        ind = np.random.choice(range(set.size(0)),bootstrap_size)
        bootstrap_set = set[ind,:]
        _, corr_mat = get_correlation_matrix(bootstrap_set, signatures, plot = False)
        corr_list.append(corr_mat)
    return np.dstack(([d for d in corr_list])).std(axis=2)

def renormalize_corr_mat(corr_mat, weights, cutoff):
    for i in range(corr_mat.shape[0]):
        for j in range(corr_mat.shape[1]):
            Ni = torch.sum(weights[:,i]>cutoff)
            Nj = torch.sum(weights[:,j]>cutoff)
            corr_mat.at[corr_mat.columns[i],corr_mat.columns[j]] = corr_mat.iloc[i,j]*np.sqrt(Ni*Nj).item()/weights.shape[0]
    return corr_mat

data_folder = "../../data/"
k_tot = 1

# Real data
real_data = csv_to_tensor("../../data/real_data/sigprofiler_not_norm_PCAWG.csv",
                              device='cpu', header=0, index_col=0)
real_data = torch.cat([real_data, torch.zeros(real_data.size(0), 7).to(real_data)], dim=1)

# Create partitions
lst_weights, lst_ctype = read_data_and_partitions(k_tot)

# Create inputs associated to the labels:
signatures = read_signatures(
    "../../data/data.xlsx", mutation_type_order="../../data/mutation_type_order.xlsx")
data_generator = DataGenerator(signatures=signatures,
                                seed=None,
                                shuffle=True)

models_path = "../../trained_models/crossval_oversample_fix"

# N_oversample_list = [1, 5, 10, 20, 30, 33]
N_oversample_list = [1]#, 10, 30]

mse_corr_oversample = []
mse_guess_oversample = []
FN_oversample = []
FP_oversample = []
for N_oversample in N_oversample_list:
    mse_corr_oversample_k = []
    mse_guess_oversample_k = []
    FN_oversample_k = []
    FP_oversample_k = []
    for k in range(k_tot):
        oversampler = CancerTypeOverSampler(lst_weights[k], lst_ctype[k])
        test_set = oversampler.get_N_oversampled_set(N_oversample)

        # Create pairs input-label
        print("Creating train, val and test data")
        test_input, test_label = data_generator.make_input(test_set, "test", "", normalize=True)

        # Run Baseline
        print("Running Baseline")
        signatures = read_signatures("../../data/data.xlsx")
        sf = Baseline(signatures)
        test_baseline = sf.get_weights_batch(test_input, n_workers=2)

        # Apply model to test set
        finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "/finetuner_low_crossval_" + str(k) + "_" + str(N_oversample),
                                        large_mum_mut_dir=models_path + "/finetuner_large_crossval_" + str(k) + "_" + str(N_oversample))
        test_guess = finetuner(mutation_dist=test_input,
                                baseline_guess=test_baseline,
                                num_mut=test_label[:, -1].view(-1, 1))

        signatures = sort_signatures("../../data/data.xlsx")
        
        # Bootstrap std
        # std_bootstrap = np.nan_to_num(bootstrap_std(test_label[:,:-1], 1000, signatures))
        std_bootstrap = np.nan_to_num(bootstrap_std(real_data, test_label.size(0), 1000, signatures))
        # pd.DataFrame(std_bootstrap).to_csv('std_bootstrap.csv')

        _, test_guess_corr = get_correlation_matrix(test_guess, signatures, plot=False)
        # _, test_label_corr = get_correlation_matrix(real_data, signatures, plot=False)
        _, test_label_corr = get_correlation_matrix(test_label[:, :-1], signatures, plot=False)
        plt.show()

        test_guess_corr = renormalize_corr_mat(test_guess_corr, test_guess, 0.01)
        test_label_corr = renormalize_corr_mat(test_label_corr, test_label[:, :-1], 0.01)
        # dif = (torch.nan_to_num(torch.Tensor(test_guess_corr.values))-torch.nan_to_num(torch.Tensor(test_label_corr.values)))**2

        # fig = plt.figure()
        # sn.heatmap(dif, annot=False)

        # fig = plt.figure()
        # sn.heatmap(std_bootstrap, annot=False)
        # plt.show()

        # std_bootstrap = torch.ones(test_guess_corr.shape)
        mse_corr_oversample_k.append(torch.sum((torch.nan_to_num(torch.Tensor(test_guess_corr.values))-torch.nan_to_num(torch.Tensor(test_label_corr.values)))**2/std_bootstrap).item())
        mse_guess_oversample_k.append(get_MSE(test_guess, test_label[:, :-1]))

        cutoff = 0.01
        label_mask = (test_label[:, :-1] > cutoff).type(torch.int).float()
        prediction_mask = (test_guess > cutoff).type(torch.int).float()
        fp = torch.sum(label_mask - prediction_mask < -0.1)
        fn = torch.sum(label_mask - prediction_mask > 0.1)
        FP_oversample_k.append(fp/test_guess.shape[0])
        FN_oversample_k.append(fn/test_guess.shape[0])

        # test_guess[test_guess<0.01] = 0
        # test_label[test_label<0.01] = 0
        # _, test_guess_corr = get_correlation_matrix(test_guess, signatures, plot=False)
        # _, test_label_corr = get_correlation_matrix(test_label[:, :-1], signatures, plot=False)
        # mse_oversample_001.append(get_MSE(torch.nan_to_num(torch.Tensor(test_guess_corr.values)), torch.nan_to_num(torch.Tensor(test_label_corr.values))).item())

        # test_guess[test_guess<0.05] = 0
        # test_label[test_label<0.05] = 0
        # _, test_guess_corr = get_correlation_matrix(test_guess, signatures, plot=False)
        # _, test_label_corr = get_correlation_matrix(test_label[:, :-1], signatures, plot=False)
        # mse_oversample_005.append(get_MSE(torch.nan_to_num(torch.Tensor(test_guess_corr.values)), torch.nan_to_num(torch.Tensor(test_label_corr.values))).item())

        # test_guess[test_guess<0.1] = 0
        # test_label[test_label<0.1] = 0
        # _, test_guess_corr = get_correlation_matrix(test_guess, signatures, plot=False)
        # _, test_label_corr = get_correlation_matrix(test_label[:, :-1], signatures, plot=False)
        # mse_oversample_01.append(get_MSE(torch.nan_to_num(torch.Tensor(test_guess_corr.values)), torch.nan_to_num(torch.Tensor(test_label_corr.values))).item())
        # plt.show()
        # plt.close()
    mse_corr_oversample.append(np.mean(mse_corr_oversample_k))
    mse_guess_oversample.append(np.mean(mse_guess_oversample_k))
    FN_oversample.append(np.mean(FN_oversample_k))
    FP_oversample.append(np.mean(FP_oversample_k))

print(mse_corr_oversample)
plt.plot(N_oversample_list, mse_corr_oversample,'--o')
plt.xlabel("Oversample")
plt.ylabel("MSE corr normalized")
plt.show()

plt.plot(N_oversample_list, mse_guess_oversample,'--o')
plt.xlabel("Oversample")
plt.ylabel("MSE guesses")
plt.show()

plt.plot(N_oversample_list, FP_oversample,'--o')
plt.xlabel("Oversample")
plt.ylabel("FP")
plt.show()

plt.plot(N_oversample_list, FN_oversample,'--o')
plt.xlabel("Oversample")
plt.ylabel("FN")
plt.show()
import itertools
import os
import sys

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import torch

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.modules.signet_module import SigNet
from signaturesnet.utilities.io import read_methods_guesses
from signaturesnet.utilities.plotting import final_plot_all_metrics_vs_mutations, final_plot_interval_metrics_vs_mutations, plot_metric_vs_mutations_classifier, plot_time_vs_mutations


def read_data():
    inputs = pd.read_csv(DATA + "/datasets/test_input.csv", header=None, index_col=None)
    labels = pd.read_csv(DATA + "/datasets/test_label.csv", header=None, index_col=None)
    num_mut = labels[72]

    signatures = pd.read_excel(DATA + '/data.xlsx')

    input_df = inputs.mul(num_mut, 0)
    input_df.index = ['sample_' + str(i) for i in list(range(len(input_df.index)))]
    input_df.columns = signatures['Type']
    print("Data loaded!")
    return inputs, input_df, labels, num_mut

if __name__ == "__main__":
    # Load data
    inputs, input_df, labels, num_mut = read_data()
    print(num_mut)

    # Load model
    signet = SigNet(classifier=TRAINED_MODELS + "/detector",
                    finetuner_realistic_low=TRAINED_MODELS + "/finetuner_low",
                    finetuner_realistic_large=TRAINED_MODELS + "/finetuner_large",
                    errorfinder=TRAINED_MODELS + "/errorfinder",
                    opportunities_name_or_path=None,
                    signatures_path=DATA + "/data.xlsx",
                    mutation_type_order=DATA + "/mutation_type_order.xlsx")
    print("model read")

    result = signet(input_df, numpy=False)
    print("forwarded")

    finetuner_guess, lower_bound, upper_bound, classification, normalized_input = result.get_output(format="numpy")
    nummuts = list(set(list(num_mut.values)))

    # guesses = finetuner_guess[:, :-1]
    # truth = labels.values[:, :-1]
    # error = (truth - guesses)
    # plt.scatter(error[:, 2], error[:, 4], s=2)
    # plt.show()
    # plt.cla()

    guesses = finetuner_guess[:, :-1]
    truth = labels.values[:, :-1]
    error = (truth - guesses) / np.sum(truth > 0.1, axis=0)
    m_corr, m_sum = np.zeros((72, 72)), np.zeros((72, 72))
    for i, j in itertools.combinations(range(72), 2):
        corr = pearsonr(error[:, i], error[:, j])
        sums = np.sum(error[:, i] + error[:, j])
        # m_corr[i, j] = corr[0]
        m_sum[i, j] = sums
    # plt.imshow(m_corr, interpolation='nearest')
    # plt.show()
    # plt.cla()
    plt.imshow(m_sum, interpolation='nearest')
    plt.show()
    print("DONE")
    # print("sig:", i, "sig:", j, ". Corr:", corr[0])


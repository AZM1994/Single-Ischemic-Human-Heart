import logging
import os
from turtle import title

import numpy as np
import pandas as pd
import torch
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.models import Generator
from signaturesnet.utilities.io import csv_to_tensor, read_model, tensor_to_csv
from signaturesnet.utilities.plotting import get_correlation_matrix, plot_correlation_matrix

device = "cuda" if torch.cuda.is_available() else "cpu"
real_data_path = os.path.join(DATA, "real_data/sigprofiler_not_norm_PCAWG.csv")
real_data_pd = pd.read_csv(real_data_path)
zero_cols = []
zero_indexes = []

def get_real_correlation_matrix():
    real_data_tensor = csv_to_tensor(real_data_path, device, header=0, index_col=0)
    return plot_correlation_matrix(real_data_tensor, sig_names=real_data_pd.columns[1:], name="real")

def get_generator_correlation_matrix():
    generator_path = os.path.join(TRAINED_MODELS, "generator")
    generator = read_model(generator_path, device)

    synt_labels = generator.generate(10000, std=1.0, cutoff=0.001)
    for i in zero_indexes:
        synt_labels[:, i] = 0
    return plot_correlation_matrix(synt_labels[:, :65], sig_names=real_data_pd.columns[1:], name="generator")

def get_synsiggen_data():
    data_file = "SynSigGen_labels.csv"
    df = pd.read_csv(data_file, index_col=0, header=0)
    # print(real_data_pd.columns[1:])
    # print(df.head)
    # print(df.describe())
    for indx, col in enumerate(real_data_pd.columns[1:]):
        if col not in df.columns:
            df.insert(indx, col, 0)
            zero_cols.append(col)
            zero_indexes.append(indx)
    print(real_data_pd.columns[1:])
    return plot_correlation_matrix(torch.tensor(df.values), sig_names=real_data_pd.columns[1:], name="synsiggen")


def get_metrics(real, guess):
    real_ = np.nan_to_num(x=real, nan=0.0, posinf=0.0, neginf=0.0)
    guess_ = np.nan_to_num(x=guess, nan=0.0, posinf=0.0, neginf=0.0)
    mask_ = np.abs(real_) > 0.05

    plt.plot((real_ - guess_)[:, 6])
    plt.show()

    diff = np.abs(real_ - guess_) * mask_

    print("diff count:", np.mean(diff > 0.01))
    print("diff mean:", np.mean(diff))
    return diff


if __name__ == "__main__":
    real = get_real_correlation_matrix()
    synsiggen = get_synsiggen_data()
    signetgen = get_generator_correlation_matrix()

    # for col in zero_cols:
    #     signetgen[col].values[:] = 0

    # signetgen_metrics = get_metrics(real, signetgen)
    # synsiggen_metrics = get_metrics(real, synsiggen)

    # plt.matshow(signetgen_metrics)
    # plt.show()

    plt.matshow(synsiggen_metrics)
    plt.show()


    sigs_names = real_data_pd.columns[1:]
    # Generator data
    generator_path = os.path.join(TRAINED_MODELS, "generator")
    generator = read_model(generator_path, device)

    synt_labels = generator.generate(1000, std=1.0)
    for i in zero_indexes:
        synt_labels[:, i] = 0
    df = pd.DataFrame(synt_labels[:, :65].cpu().detach().numpy(), columns=sigs_names)
    corrMatrix_gen = df.corr()

    # SynSigGen data
    data_file = "SynSigGen_labels.csv"
    df = pd.read_csv(data_file, index_col=0, header=0)
    for indx, col in enumerate(real_data_pd.columns[1:]):
        if col not in df.columns:
            df.insert(indx, col, 0)
            zero_cols.append(col)
            zero_indexes.append(indx)
    df = pd.DataFrame(df, columns=sigs_names)
    corrMatrix_synsiggen = df.corr()

    # Real data
    real_data_tensor = csv_to_tensor(real_data_path, device, header=0, index_col=0)
    df = pd.DataFrame(real_data_tensor.cpu().detach().numpy(), columns=sigs_names)
    corrMatrix_real = df.corr()

    # Plot
    import seaborn as sn
    fig,ax = plt.subplots(1,4, figsize=(16,5), gridspec_kw=dict(width_ratios=[1,1,1,0.07]))
    sn.heatmap(corrMatrix_real, annot=False, vmin=-0.6, vmax=1, cbar=False, ax=ax[0])
    ax[0].set_title('Real correlations')
    sn.heatmap(corrMatrix_synsiggen, annot=False, vmin=-0.6, vmax=1, cbar=False, ax=ax[1])
    ax[1].set_title('SynSigGen correlations')
    sn.heatmap(corrMatrix_gen, annot=False, vmin=-0.6, vmax=1, cbar=False, ax=ax[2])
    ax[2].set_title('SigNet Generator correlations')
    fig.colorbar(ax[2].collections[0], cax=ax[3])
    fig.tight_layout()
    plt.savefig('generator_correlations.pdf')


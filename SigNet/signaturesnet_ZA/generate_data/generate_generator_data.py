import torch
import matplotlib.pyplot as plt
import numpy as np

from signaturesnet.utilities.io import read_model, read_data_generator

def make_hist(data):
    indexes = []
    for example in data:
        idx = (example > 0).nonzero().flatten()
        indexes += idx.tolist()
    plt.hist(indexes, bins=72)
    plt.show()

def make_hist_comparison(real, data):
    indexes = []
    for example in data:
        idx = (example > 0).nonzero().flatten()
        indexes += idx.tolist()
    indexes_real = []
    for example in real:
        idx = (example > 0).nonzero().flatten()
        indexes_real += idx.tolist()

    counts_real, bins, _ = plt.hist(indexes_real, bins=range(72), density=True, align="left", rwidth=0.4)
    counts, _, _ = plt.hist(indexes, bins=bins, density=True, align="mid", rwidth=0.4)
    plt.legend(['Real', 'Generated'])
    plt.ylabel('Density')
    plt.show()

def make_mean_weight_comparison(real, data):
    mean_weights_real = torch.mean(real, dim=0).detach().numpy()
    mean_weights = torch.mean(data, dim=0).detach().numpy()

    fig, ax = plt.subplots()
    ax.bar(range(72),mean_weights_real, align='center', width=0.2, alpha=0.5, ecolor='black', capsize=10)
    ax.bar(np.array(range(72))+0.2, mean_weights, width=0.2, align='center')
    ax.set_ylabel('Mean weight')
    ax.set_xlabel('Signature')
    plt.legend(['Real', 'Generated'])
    plt.show()

if __name__ == "__main__":
    generator = read_model("../trained_models/exp_good/generator", device="cpu")

    examples = generator.generate(5000)
    examples[examples <= 0.01] = 0
    # make_hist(examples)

    real_data, val_real_data = read_data_generator(device="cpu")
    real_data.inputs[real_data.inputs <= 0.01] = 0
    # make_hist(real_data.inputs)

    make_hist_comparison(real_data.inputs, examples)
    make_mean_weight_comparison(real_data.inputs, examples)
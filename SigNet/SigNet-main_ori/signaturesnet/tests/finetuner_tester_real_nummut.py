import numpy as np
import pandas as pd
from signaturesnet.utilities.metrics import get_classification_metrics
import torch

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.models import Baseline
from signaturesnet.modules import CombinedFinetuner
from signaturesnet.utilities.io import csv_to_tensor, read_signatures
from signaturesnet.utilities.plotting import plot_all_metrics_vs_mutations

def sample_from_sig(signature, num_mut):
        sample = signature
        c = torch.distributions.categorical.Categorical(probs=signature)
        sample_shape = torch.Size([num_mut, ])
        samples = c.sample(sample_shape=sample_shape).type(torch.float32)
        n_bins = signature.shape[0]
        sample = torch.histc(samples, bins=n_bins, min=0, max=n_bins - 1)
        sample = sample/float(num_mut)
        return sample

if __name__ == "__main__":
    real_data_weights = pd.read_csv(DATA + "/real_data/sigprofiler_not_norm_PCAWG.csv", header=0, index_col=0).reset_index(drop=True)
    real_data_weights = pd.concat([real_data_weights, pd.DataFrame(np.zeros((real_data_weights.shape[0], 7)))], axis=1, ignore_index=True)
    real_data_weights = torch.tensor(real_data_weights.values, dtype=torch.float)

    real_data = pd.read_csv(DATA + '/real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv', sep=',')
    real_data = real_data.iloc[:,3:]
    total_muts = torch.tensor(real_data.sum(axis=1), dtype=torch.float)

    signatures = read_signatures("../data/data.xlsx", mutation_type_order="../data/mutation_type_order.xlsx")

    n_samples = 100
    inputs =  torch.empty((real_data.shape[0]*n_samples, 96))
    for i, num_mut in enumerate(total_muts):
        # Compute resulting signature
        signature = torch.einsum("ij,j->i", (signatures, real_data_weights[i]))
        for j in range(n_samples):
            # Sample the mutational vector from the given distribution
            sample = sample_from_sig(signature=signature,
                                        num_mut=int(num_mut))
            # Store Values
            inputs[i*n_samples+j, :] = sample
    total_muts = torch.tensor([n for n in total_muts for _ in range(n_samples)], dtype=torch.float)
    total_real_data_weights = real_data_weights.repeat_interleave(n_samples, 0)
    # Run Baseline
    print("Running Baseline")
    sf = Baseline(signatures)
    # test_baseline =  torch.empty((inputs.shape[0], 72))
    test_baseline = sf.get_weights_batch(inputs)
    finetuner = CombinedFinetuner(low_mum_mut_dir=TRAINED_MODELS + "/finetuner_low",
                                    large_mum_mut_dir=TRAINED_MODELS + "/finetuner_large")
    finetuner_guess = finetuner(mutation_dist=inputs,
                                baseline_guess=test_baseline,
                                num_mut=total_muts,
                                cutoff_0=0.01)
    input_batch = csv_to_tensor(DATA + "/datasets/test_input.csv")
    label_batch = csv_to_tensor(DATA + "/datasets/test_label.csv")

    # Run Baseline
    print("Running Baseline")
    test_baseline = sf.get_weights_batch(input_batch, n_workers=4)
    finetuner_test = finetuner(mutation_dist=input_batch,
                                baseline_guess=test_baseline,
                                num_mut=label_batch[:, -1].view(-1, 1),
                                cutoff_0=0.01)

    import matplotlib.pyplot as plt
    # metrics = get_classification_metrics(label_batch=real_data_weights,
    #                                     prediction_batch=finetuner_guess[:, :-1])
    # plt.plot(total_muts, metrics['accuracy %'])

    num_muts = np.unique(label_batch[:,-1].detach().numpy())
    MAE_all_test = [0]*10
    MAE_real = [0]*10
    for i, num_mut in enumerate(num_muts):
        indexes = label_batch[:, -1] == num_mut
        MAE_all_test[i] = torch.mean(torch.abs(label_batch[indexes, :-1]-finetuner_test[indexes, :-1])).item()
        if i < len(num_muts)-1 and i > 0:
            indexes = (total_muts > (num_mut-num_muts[i-1])/2) & (total_muts < (num_muts[i+1]-num_mut)/2) 
            print(torch.sum(indexes).item())
        if i == 0:
            indexes = (total_muts < (num_muts[i+1]-num_mut)/2)
            print(torch.sum(indexes).item())
        if i == len(num_muts)-1:
            indexes = (total_muts > (num_mut-num_muts[i-1])/2)
            print(torch.sum(indexes).item())

        MAE_real[i] = torch.mean(torch.abs(total_real_data_weights[indexes, :]-finetuner_guess[indexes, :-1])).item()
    print(MAE_all_test)
    print(MAE_real)
    # plt.bar(np.log10(num_muts), MAE_all_test, alpha=0.4)
    # plt.scatter(np.log10(total_muts), torch.mean(torch.abs(real_data_weights - finetuner_guess[:, :-1]),dim=1), alpha=0.5)
    # plt.show()

    plt.scatter(np.log10(num_muts), MAE_all_test, alpha=0.4)
    plt.scatter(np.log10(num_muts), MAE_real, alpha=0.5)
    plt.legend(['Test set', 'Real set'])
    plt.ylabel('Mean absolute error')
    plt.xlabel('log(N)')
    plt.show()
    

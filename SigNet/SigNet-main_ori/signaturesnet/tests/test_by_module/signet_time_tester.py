import os

import torch
import time
import numpy as np
import pandas as pd

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.models.baseline import Baseline
from signaturesnet.modules.combined_finetuner import CombinedFinetuner
from signaturesnet.utilities.io import csv_to_tensor, read_model, read_signatures

# Read data

# Load data
data_path = "../../../data/exp_all/"
inputs = csv_to_tensor(data_path + "test_input.csv", device='cpu')
labels = csv_to_tensor(data_path + "test_label.csv", device='cpu')
num_mut = labels[:, -1].unique()
print("data loaded")

# Paths
signatures_path = os.path.join(DATA, "data.xlsx")
mutation_type_order = os.path.join(DATA, "mutation_type_order.xlsx")
finetuner_realistic_low = os.path.join(TRAINED_MODELS, "finetuner_low")
finetuner_realistic_large = os.path.join(TRAINED_MODELS, "finetuner_large")
classifier = os.path.join(TRAINED_MODELS, "detector")
errorfinder = os.path.join(TRAINED_MODELS, "errorfinder")

# Load models
signatures = read_signatures(file=signatures_path,
                             mutation_type_order=mutation_type_order)
baseline = Baseline(signatures)
finetuner = CombinedFinetuner(low_mum_mut_dir=finetuner_realistic_low,
                                      large_mum_mut_dir=finetuner_realistic_large)
classifier = read_model(classifier)
errorfinder = read_model(errorfinder)

list_of_inputs = []
for mut in num_mut:
    list_of_inputs.append(inputs[labels[:, -1] == mut, :])

# Functions
def join_and_sort(low, large, ind_order):
        joined = torch.cat((low, large), dim=0)
        joined = torch.cat((joined, ind_order), dim=1)
        joined = joined[joined[:, -1].sort()[1]]
        return joined[:, :-1]

def apply_cutoff(comb, cutoff):
        mask = (comb > cutoff).type(torch.int).float()
        comb = comb*mask
        comb = torch.cat((comb, torch.ones_like(torch.sum(
        comb, axis=1).reshape((-1, 1)))-torch.sum(
        comb, axis=1).reshape((-1, 1))), axis=1)
        return comb

def separate_classification(classification, mutation_dist, baseline_guess, num_mut):
    classification_cutoff = 0.5
    ind = torch.tensor(range(classification.size()[0]))
    ind_order = torch.tensor(np.concatenate((ind[classification <= classification_cutoff], ind[classification > classification_cutoff]))).reshape(-1, 1).to(torch.float).to('cpu')
    
    input_batch_random = mutation_dist[classification <= classification_cutoff, ]
    input_batch_realistic = mutation_dist[classification > classification_cutoff, ]
    num_mut_realistic = num_mut[classification > classification_cutoff, ]
    classification_realistic = classification[classification > classification_cutoff, ]

    baseline_guess_random = baseline_guess[classification <= classification_cutoff, ]
    baseline_guess_realistic = baseline_guess[classification > classification_cutoff, ]

    baseline_guess_random = baseline_guess_random/torch.sum(baseline_guess_random, dim=1).reshape(-1,1)
    
    return input_batch_realistic, input_batch_random, baseline_guess_random, baseline_guess_realistic, num_mut_realistic, classification_realistic, ind_order


replicates = 5
times_baseline = np.zeros((replicates, len(num_mut)))
times_classifier = np.zeros((replicates, len(num_mut)))
times_finetuner = np.zeros((replicates, len(num_mut)))
times_errorfinder = np.zeros((replicates, len(num_mut)))
times_total = np.zeros((replicates, len(num_mut)))
for i,input in enumerate(list_of_inputs):
    for k in range(replicates):
        st_total = time.time()
        st = time.time()
        baseline_guess = baseline.get_weights_batch(input_batch=input, 
                                                    n_workers=1)
        et = time.time()
        times_baseline[k, i] = (et-st)/baseline_guess.size(dim=0)                                                                
        
        st = time.time()
        num_muts_i = num_mut[i].repeat(input.size(dim=0),1)
        classification = classifier(mutation_dist=input,
                                    num_mut=num_muts_i).view(-1)
        et = time.time()
        times_classifier[k, i] = (et-st)/classification.size(dim=0)     

        mutation_dist_realistic, mutation_dist_random, baseline_guess_random, baseline_guess_realistic, num_mut_realistic, classification_realistic, ind_order = separate_classification(classification, input, baseline_guess, num_muts_i)

        st = time.time()
        finetuner_guess_realistic = finetuner(mutation_dist=mutation_dist_realistic,
                                         baseline_guess = baseline_guess_realistic,
                                         num_mut=num_mut_realistic)

        baseline_guess_random = apply_cutoff(baseline_guess_random, 0.01)
        finetuner_guess = join_and_sort(finetuner_guess_realistic, baseline_guess_random, ind_order)
        et = time.time()
        times_finetuner[k, i] = (et-st)/finetuner_guess.size(dim=0) 

        st = time.time()
        upper, lower = errorfinder(weights=finetuner_guess_realistic[:,:-1],
                                           num_mutations=num_mut_realistic,
                                           classification=classification_realistic.reshape(-1, 1))

        upper = join_and_sort(upper, torch.full_like(baseline_guess_random[:,:-1], float('nan')), ind_order)
        lower = join_and_sort(lower, torch.full_like(baseline_guess_random[:,:-1], float('nan')), ind_order)
        et = time.time()
        times_errorfinder[k, i] = (et-st)/lower.size(dim=0)  

        et_total = time.time()
        times_total[k, i] = (et_total-st_total)/input.size(dim=0) 


times_df = pd.DataFrame(times_baseline)
times_df.columns = num_mut.tolist()
times_df.loc['mean'] = times_df.mean()
times_df.to_csv('baseline_times.txt', index=False)

times_df = pd.DataFrame(times_classifier)
times_df.columns = num_mut.tolist()
times_df.loc['mean'] = times_df.mean()
times_df.to_csv('classifier_times.txt', index=False)

times_df = pd.DataFrame(times_finetuner)
times_df.columns = num_mut.tolist()
times_df.loc['mean'] = times_df.mean()
times_df.to_csv('finetuner_times.txt', index=False)

times_df = pd.DataFrame(times_errorfinder)
times_df.columns = num_mut.tolist()
times_df.loc['mean'] = times_df.mean()
times_df.to_csv('errorfinder_times.txt', index=False)

times_df = pd.DataFrame(times_total)
times_df.columns = num_mut.tolist()
times_df.loc['mean'] = times_df.mean()
times_df.to_csv('total_times.txt', index=False)

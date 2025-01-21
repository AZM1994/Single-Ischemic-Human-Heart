import os
import sys

import pandas as pd


sys.path.append(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
from models import Baseline, ErrorFinder
from modules import CombinedFinetuner
from utilities.io import csv_to_tensor, read_model, read_signatures
from utilities.plotting import plot_all_metrics_vs_mutations,\
                               plot_metric_vs_mutations,\
                               plot_interval_metrics_vs_mutations,\
                               plot_interval_performance

# Load data
data_path = "../../../data/exp_generator/"
inputs = csv_to_tensor(data_path + "test_generator/test_generator_input.csv", device='cpu')
# baselines = csv_to_tensor(data_path + "test_generator/test_generator_baseline.csv", device='cpu')
labels = csv_to_tensor(data_path + "test_generator/test_generator_label.csv", device='cpu')
num_mut = labels[:, -1].unsqueeze(1)

signatures = read_signatures(file="../../../data/data.xlsx",
                             mutation_type_order="../../../data/mutation_type_order.xlsx")
baseline = Baseline(signatures)
baselines = baseline.get_weights_batch(inputs)

print(baselines)

# Make guess
models_path = "../../../trained_models/exp_generator/"
finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "finetuner_generator_low_2",
                              large_mum_mut_dir=models_path + "finetuner_generator_large")
finetuner_guess = finetuner(mutation_dist=inputs,
                            baseline_guess=baselines,
                            num_mut=num_mut)

# Plot results
list_of_methods = ['baseline', 'finetuner']
list_of_guesses = [baselines, finetuner_guess]
# plot_all_metrics_vs_mutations( list_of_methods, list_of_guesses, labels, '')
plot_metric_vs_mutations(list_of_metrics=["accuracy %", "reconstruction_error"],
                         list_of_methods=list_of_methods,
                         list_of_guesses=list_of_guesses,
                         label=labels,
                         show=True,
                         signatures=signatures,
                         mutation_distributions=inputs)

classifier = read_model('../../../trained_models/exp_generator/classifier')
classifier_output = classifier(mutation_dist=inputs,
                               num_mut=num_mut)
# # Run errorfinder
errorfinder = read_model(models_path + "errorfinder_generator_103")
upper_guess, lower_guess = errorfinder(weights=finetuner_guess,
                                       num_mutations=num_mut,
                                       classification=classifier_output)

# # Plot results
plot_interval_metrics_vs_mutations(labels, upper_guess, lower_guess, '', True)
plot_interval_performance(labels, upper_guess, lower_guess, list(pd.read_excel("../../../data/data.xlsx").columns)[1:], '', True)

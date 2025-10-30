import os
import sys

import pandas as pd
import torch
import numpy as np

from utilities.plotting import plot_metric_vs_mutations, plot_metric_vs_sigs

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.io import read_methods_guesses
from utilities.metrics import get_classification_metrics
from models import Baseline, FineTuner

dev = torch.device("cpu")

experiment_id = "exp_0"
# Our methods (baseline with YAPSA and finetuner)
num_classes = 72
signatures_data = pd.read_excel("../data/data.xlsx")
signatures = [torch.tensor(signatures_data.iloc[:, i]).type(torch.float32)
                  for i in range(2, 74)][:num_classes]
input_batch = torch.tensor(pd.read_csv(
        "../data/%s/test_realistic/test_realistic_input.csv"%experiment_id, header=None).values, dtype=torch.float)
label_batch = torch.tensor(pd.read_csv(
        "../data/%s/test_realistic/test_realistic_label.csv"%experiment_id, header=None).values, dtype=torch.float)

sf = Baseline(signatures)
baseline_guess = sf.get_weights_batch(input_batch)

finetuner_model_name = "finetuner_mixed"
finetuner_params = {"num_hidden_layers": 3,
                        "num_units": 600,
                        "num_classes": 72}

finetuner = FineTuner(**finetuner_params)
finetuner.load_state_dict(torch.load(os.path.join(
            "../trained_models/%s"%experiment_id, finetuner_model_name), map_location=torch.device('cpu')))
finetuner.eval() # NOTE: Important! Otherwise we don't zero small values
finetuner_guess = finetuner(mutation_dist=input_batch, weights=baseline_guess, num_mut=label_batch[:,-1].reshape(-1,1))


list_of_methods = ["decompTumor2Sig", "MutationalPatterns", "mutSignatures", "SignatureEstimationQP","YAPSA"] # "deconstructSigs"

list_of_guesses, label = read_methods_guesses(dev, experiment_id, "test_realistic", list_of_methods)

list_of_guesses.append(baseline_guess)
list_of_guesses.append(finetuner_guess)
list_of_methods.append("Baseline")
list_of_methods.append("Finetuner")

list_of_metrics = ["MAE_p", "MAE_n", "fp", "fn"]

plot_metric_vs_sigs(list_of_metrics, list_of_methods, list_of_guesses, label, "mixed_realistic_vs_sigs")
plot_metric_vs_mutations(list_of_metrics, list_of_methods, list_of_guesses, label, "mixed_realistic")

list_of_metrics = ["accuracy %", "sens: tp/p %", "spec: tn/n %"]

plot_metric_vs_sigs(list_of_metrics, list_of_methods, list_of_guesses, label, "mixed_realistic_performance_vs_sigs")
plot_metric_vs_mutations(list_of_metrics, list_of_methods, list_of_guesses, label, "mixed_realistic_performance")


# for method in list_of_methods:
#     print(method)
#     metrics = get_classification_metrics(label_batch=label[..., :-1], prediction_batch=globals()[method + "_guess_random"])
#     for metric in metrics:
#         print(metric, np.round(metrics[metric].detach().numpy(), 8))
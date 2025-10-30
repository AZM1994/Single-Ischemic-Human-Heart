import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import seaborn as sn


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from models import Baseline, ErrorFinder, FineTuner
from utilities.io import read_test_data
from utilities.metrics import *
from utilities.plotting import plot_interval_metrics_vs_mutations, plot_interval_metrics_vs_sigs, plot_interval_width_vs_mutations, plot_interval_width_vs_mutations_some_sigs, plot_confusion_matrix, plot_weights, plot_weights_comparison, plot_interval_performance


class ModelTester:
    def __init__(self, num_classes):
        self.num_classes = num_classes

    def get_basic_metrics(self, guessed_labels, true_labels):
            metrics = {
                "mse": get_MSE,
                "cos": get_negative_cosine_similarity,
                "KL": get_kl_divergence,
                "JS": get_jensen_shannon,
                "W": get_wasserstein_distance,
            }
            for key in metrics.keys():
                val = metrics[key](predicted_label=guessed_labels,
                                true_label=true_labels)
                print(key, ":", np.round(val, decimals=3))

    def get_confusion_matrix(self, guessed_labels, true_labels):
        label_sigs, predicted_sigs = self.__probs_batch_to_sigs(
            true_labels, guessed_labels, cutoff=0.05, num_classes=self.num_classes)

        conf_mat = plot_confusion_matrix(
            label_sigs, predicted_sigs, range(self.num_classes+1))

    def __probs_batch_to_sigs(self, label_batch, predicted_batch, cutoff=0.05, num_classes=72):
        label_sigs_list = list(range(num_classes))
        predicted_sigs_list = list(range(num_classes))
        for i in range(len(label_batch)):
            for j in range(len(label_batch[i])):
                if label_batch[i][j] > cutoff and predicted_batch[i][j] > cutoff:
                    label_sigs_list.append(j)
                    predicted_sigs_list.append(j)
                    continue
                if label_batch[i][j] > cutoff and predicted_batch[i][j] < cutoff:
                    label_sigs_list.append(j)
                    predicted_sigs_list.append(num_classes)
                    continue
                if label_batch[i][j] < cutoff and predicted_batch[i][j] > cutoff:
                    label_sigs_list.append(num_classes)
                    predicted_sigs_list.append(j)
                    continue
        return torch.tensor(label_sigs_list), torch.tensor(predicted_sigs_list)


if __name__ == "__main__":
    experiment_id = "exp_0"
    test_id = "test_realistic"
    device = "cpu"

    # Model params finetuner
    model_id_finetuner = "finetuner_mixed_js_loss"
    num_hidden_layers = 2
    num_neurons = 600
    num_classes = 72

    # Model params error
    model_id_error_learner = "errorfinder_mixed"
    num_hidden_layers_pos = 2
    num_neurons_pos = 1000
    num_hidden_layers_neg = 2
    num_neurons_neg = 1000

    # Open data
    data = pd.read_excel("../../data/data.xlsx")
    signatures = [torch.tensor(data.iloc[:, i]).type(torch.float32)
                  for i in range(2, 74)][:num_classes]

    input_batch, label_mut_batch = read_test_data(device=device,
                                                  experiment_id=experiment_id,
                                                  test_id=test_id,
                                                  data_folder="../../data")
    n_datapoints = -1
    input_batch = input_batch[:n_datapoints]
    label_batch = label_mut_batch[:n_datapoints, :-1]
    num_mut = label_mut_batch[:n_datapoints, -1].reshape((-1, 1))

    # Baseline:
    sf = Baseline(signatures)
    baseline_batch = sf.get_weights_batch(input_batch)  # [:50, ...])

    # Instantiate model and do predictions for finetuner:
    model = FineTuner(num_classes=num_classes,
                      num_hidden_layers=num_hidden_layers,
                      num_units=num_neurons)
    model.load_state_dict(torch.load(os.path.join(
        "../../trained_models/" + experiment_id, model_id_finetuner), map_location=torch.device(device)))
    model.eval()
    guessed_labels = model(input_batch, baseline_batch, num_mut)

    # Instantiate model and do predictions for error learner:
    model_error = ErrorFinder(num_classes=num_classes,
                              num_hidden_layers_pos=num_hidden_layers_pos,
                              num_units_pos=num_neurons_pos,
                              num_hidden_layers_neg=num_hidden_layers_neg,
                              num_units_neg=num_neurons_neg)

    model_error.load_state_dict(torch.load(os.path.join(
        "../../trained_models/" + experiment_id, model_id_error_learner), map_location=torch.device(device)))
    model_error.eval()
    pred_upper, pred_lower = model_error(guessed_labels, num_mut)

    # # Get metrics
    # model_tester = ModelTester(num_classes=num_classes)
    # model_tester.get_basic_metrics(guessed_labels=guessed_labels, true_labels=label_batch)

    # # Plot signatures
    plot_weights_comparison(label_batch[0, :].detach().numpy(), guessed_labels[0, :].detach().numpy(
    ), pred_upper[0, :].detach().numpy(),pred_lower[0, :].detach().numpy(), list(data.columns)[2:], "example_25mut.png")
    plot_weights_comparison(label_batch[5000, :].detach().numpy(), guessed_labels[5000, :].detach().numpy(
    ), pred_upper[5000, :].detach().numpy(),pred_lower[5000, :].detach().numpy(), list(data.columns)[2:],"example_150mut.png")
    plot_weights_comparison(label_batch[12000, :].detach().numpy(), guessed_labels[12000, :].detach().numpy(
    ), pred_upper[12000, :].detach().numpy(),pred_lower[12000, :].detach().numpy(), list(data.columns)[2:],"example_10kmut.png")


    # Plot interval performance
    plot_interval_metrics_vs_sigs(label_mut_batch, pred_upper, pred_lower, "mixed_realistic_interval_vs_sigs")
    plot_interval_metrics_vs_mutations(label_mut_batch, pred_upper, pred_lower, "mixed_realistic_interval")
    plot_interval_performance(label_batch, pred_upper,pred_lower, list(data.columns)[2:], "mixed_realistic_interval_performance")

    # Plot interval width vs number of mutations
    plot_interval_width_vs_mutations(label_mut_batch, 
        pred_upper, pred_lower, list(data.columns)[2:])
    
    list_of_sigs_ind = [0,4,5,44]
    plot_interval_width_vs_mutations_some_sigs(label_mut_batch, pred_upper, pred_lower, list_of_sigs_ind, list(data.columns)[2:])

    # Print metrics
    # in_prop, mean_length = get_pi_metrics(
    #     label_batch, pred_lower=pred_lower, pred_upper=pred_upper)
    # print("In proportion:", in_prop.detach().numpy())
    # print("Mean length:", mean_length.detach().numpy())

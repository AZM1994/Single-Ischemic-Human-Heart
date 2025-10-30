import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from models import Baseline
from utilities.io import read_signatures, read_cosmic_v2_signatures, read_test_data, read_model
from utilities.plotting import plot_metric_vs_mutations, plot_metric_vs_sigs

data_id = "exp_v2"
test_id = "test_random"
model_id = "../../trained_models/finetuner_cosmicv2"

input_batch, label_batch = read_test_data("cpu", data_id, test_id, data_folder="../../data")
signatures = read_cosmic_v2_signatures("../../data/data_v2.xlsx")

baseline = Baseline(signatures)
baseline_guess = baseline.get_weights_batch(input_batch, n_workers=2)
label = label_batch

finetuner = read_model(directory=model_id)
finetuner_guess = finetuner(mutation_dist=input_batch,
                            baseline_guess=baseline_guess,
                            num_mut=label_batch[:,-1].view(-1, 1))


res_dict = {'baseline':baseline_guess, 'finetuner':finetuner_guess}

list_of_metrics = ["MAE_p", "MAE_n", "fp", "fn"]

plot_metric_vs_sigs(list_of_metrics, list(res_dict.keys()), list(res_dict.values()), label, "../../plots/%s/random_vs_sigs.png"%data_id)
plot_metric_vs_mutations(list_of_metrics, list(res_dict.keys()), list(res_dict.values()), label, "../../plots/%s/random_vs_muts.png"%data_id)

list_of_metrics = ["accuracy %", "sens: tp/p %", "spec: tn/n %"]

plot_metric_vs_sigs(list_of_metrics, list(res_dict.keys()), list(res_dict.values()), label, "../../plots/%s/random_accuracy_vs_sigs.png"%data_id)
plot_metric_vs_mutations(list_of_metrics, list(res_dict.keys()), list(res_dict.values()), label, "../../plots/%s/random_accuracy_vs_muts.png"%data_id)

import matplotlib.pyplot as plt
import numpy as np

x = np.log10(label[:,-1].detach().numpy())
plt.hist(x)
plt.show()

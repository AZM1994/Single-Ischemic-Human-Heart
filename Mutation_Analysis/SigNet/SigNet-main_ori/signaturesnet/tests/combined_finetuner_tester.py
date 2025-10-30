import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules import CombinedFinetuner
from models import Baseline
from utilities.io import read_methods_guesses, read_signatures, read_test_data
from utilities.plotting import plot_all_metrics_vs_mutations, plot_metric_vs_mutations, plot_metric_vs_sigs

experiment_id = "exp_v2"

test_id = "test"

large_mum_mut_dir = "../../trained_models/%s/finetuner_realistic_large_v2_03"%experiment_id
low_mum_mut_dir = "../../trained_models/%s/finetuner_realistic_low_v2_03"%experiment_id

input_batch, label_batch = read_test_data("cpu", experiment_id, test_id, data_folder="../../data")
signatures = read_signatures("../../data/data_v2.xlsx")

baseline = Baseline(signatures)
baseline_guess = baseline.get_weights_batch(input_batch, n_workers=1)

finetuner = CombinedFinetuner(low_mum_mut_dir=low_mum_mut_dir, large_mum_mut_dir=large_mum_mut_dir)

finetuner_guess = finetuner(input_batch, baseline_guess, label_batch[:,-1])

list_of_methods = ["decompTumor2Sig", "mutSignatures", "SignatureEstimationQP","YAPSA"]#, "deconstructSigs"]#, "MutationalPatterns",
list_of_guesses, label = read_methods_guesses('cpu', experiment_id, "test", list_of_methods, data_folder="../../data")

list_of_methods += ['Baseline', 'Finetuner']
list_of_guesses += [baseline_guess, finetuner_guess]

plot_all_metrics_vs_mutations(list_of_methods, list_of_guesses, label_batch, '../../plots/exp_final')

# list_of_metrics = ["MAE_p", "MAE_n", "fp", "fn"]

# plot_metric_vs_sigs(list_of_metrics, list_of_methods, list_of_guesses, label, "../../plots/%s/random_vs_sigs.png"%experiment_id)
# plot_metric_vs_mutations(list_of_metrics, list_of_methods, list_of_guesses, label, "../../plots/%s/random_vs_muts.png"%experiment_id)

# list_of_metrics = ["accuracy %", "sens: tp/p %", "spec: tn/n %"]

# plot_metric_vs_sigs(list_of_metrics, list_of_methods, list_of_guesses, label, "../../plots/%s/random_accuracy_vs_sigs.png"%experiment_id)
# plot_metric_vs_mutations(list_of_metrics, list_of_methods, list_of_guesses, label, "../../plots/%s/random_accuracy_vs_muts.png"%experiment_id)

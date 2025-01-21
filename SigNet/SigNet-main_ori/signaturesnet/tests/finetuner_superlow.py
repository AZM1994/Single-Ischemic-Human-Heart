
from signaturesnet.modules import CombinedFinetuner
from signaturesnet.utilities.io import csv_to_tensor
from signaturesnet.utilities.plotting import plot_all_metrics_vs_mutations_superlow

path = "../../data/exp_superlow_nummut/refitter/"
input_batch = csv_to_tensor(path + "val_superlow_input.csv")
label_batch = csv_to_tensor(path + "val_superlow_label.csv")
test_baseline = csv_to_tensor(path + "val_superlow_baseline.csv")

models_path = "../trained_models/"
finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "finetuner_low",
                                large_mum_mut_dir=models_path + "finetuner_large")
finetuner_guess = finetuner(mutation_dist=input_batch,
                            baseline_guess=test_baseline,
                            num_mut=label_batch[:, -1].view(-1, 1),
                            cutoff_0=0.01)

list_of_methods = ['NNLS', 'Finetuner']
list_of_guesses = [test_baseline, finetuner_guess[:,:-1]]

plot_all_metrics_vs_mutations_superlow(list_of_methods, list_of_guesses, label_batch, '', show=True)


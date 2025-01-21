import numpy as np

from signaturesnet import DATA
from signaturesnet.models import Baseline
from signaturesnet.modules import CombinedFinetuner
from signaturesnet.trainers.finetuner_crossvalidation import read_data_and_partitions
from signaturesnet.utilities.io import csv_to_tensor, read_methods_guesses, read_signatures
from signaturesnet.utilities.plotting import plot_all_metrics_vs_mutations, plot_crossval, plot_crossval_benchmark
from signaturesnet.utilities.metrics import get_classification_metrics
from signaturesnet.utilities.data_generator import DataGenerator

experiment_id = "crossval"
crossval = True

if crossval == True:

    # Plots crossvalidation
    k_tot = 10

    # Create partitions
    lst_weights, lst_ctype = read_data_and_partitions(k_tot)

    # Create inputs associated to the labels:
    signatures = read_signatures(
        DATA +"/data.xlsx", mutation_type_order=DATA +"/mutation_type_order.xlsx")
    data_generator = DataGenerator(signatures=signatures,
                                    seed=None,
                                    shuffle=True)

    models_path = "../../trained_models/" + experiment_id

    # Loop through the partitions
    list_of_metrics = ["MAE", "KL", "fpr", "fnr", "accuracy %",
                        "precision %", "sens: tp/p %", "spec: tn/n %"]
    num_muts = [25, 50, 100, 250, 500, 1e3, 5e3, 1e4, 5e4, 1e5]
    values = np.zeros((k_tot, len(num_muts), len(list_of_metrics)))
    for k in range(k_tot):
        current_test = k

        # Create test weight sets
        test_weights = lst_weights[current_test]

        # Create pairs input-label
        test_input, test_label = data_generator.make_input(test_weights, "test", "large")
        
        for i, num_mut in enumerate(num_muts):
            indexes = test_label[:, -1] == num_mut

            # Run Baseline
            print("Running Baseline")
            sf = Baseline(signatures)
            test_baseline = sf.get_weights_batch(test_input[indexes,:], n_workers=4)

            # Apply model to test set
            finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "/finetuner_low_crossval_" + str(i),
                                            large_mum_mut_dir=models_path + "/finetuner_large_crossval_" + str(i))
            finetuner_guess = finetuner(mutation_dist=test_input[indexes,:],
                                        baseline_guess=test_baseline,
                                        num_mut=test_label[indexes, -1].view(-1, 1))
            
            # Test model
            metrics = get_classification_metrics(label_batch=test_label[indexes, :-1],
                                                prediction_batch=finetuner_guess[:,:-1])
            for metric_index, metric in enumerate(list_of_metrics):
                values[k, i, metric_index] = metrics[metric]

    # Plot final results
    plot_crossval(values, num_muts)

    # Compute final metrics crossvalidation and benchmark:
    test_input = csv_to_tensor("../../data/exp_all/test_input.csv")
    # Load Baseline and get guess
    baseline = Baseline(signatures)
    baseline_guess = baseline.get_weights_batch(test_input, n_workers=2)

    list_of_methods = ["decompTumor2Sig", "MutationalPatterns", "mutSignatures", "SignatureEstimationQP","YAPSA", "deconstructSigs", "mutationalCone", "QPsig"]
    list_of_guesses, label = read_methods_guesses('cpu', "exp_all", list_of_methods, data_folder="../../data/")
    list_of_methods += ['NNLS']
    list_of_guesses += [baseline_guess]

    plot_crossval_benchmark(list_of_methods, list_of_guesses, label, values, folder_path='../../plots', show=False)


# Without crossval:
if crossval == False:

    model_names = ["_fp_1", "_fp_01", "_fp_001", ""]

    for model_name in model_names:
        input_batch = csv_to_tensor("../../data/exp_all/test_input.csv")
        label_batch = csv_to_tensor("../../data/exp_all/test_label.csv")

        signatures = read_signatures("../../data/data.xlsx")
        # Run Baseline
        print("Running Baseline")
        sf = Baseline(signatures)
        test_baseline = sf.get_weights_batch(input_batch, n_workers=4)

        models_path = "../../trained_models/exp_fp"
        finetuner = CombinedFinetuner(low_mum_mut_dir=models_path + "/finetuner_low_1" + model_name,
                                        large_mum_mut_dir=models_path + "/finetuner_large_1" + model_name)
        finetuner_guess = finetuner(mutation_dist=input_batch,
                                    baseline_guess=test_baseline,
                                    num_mut=label_batch[:, -1].view(-1, 1))

        list_of_methods = ["decompTumor2Sig", "MutationalPatterns", "mutSignatures", "SignatureEstimationQP","YAPSA"]#, "deconstructSigs"]
        list_of_guesses, label = read_methods_guesses('cpu', "exp_all", list_of_methods, data_folder="../../data/")
        list_of_methods += ['NNLS', 'Finetuner']
        list_of_guesses += [test_baseline, finetuner_guess]

        plot_all_metrics_vs_mutations(list_of_methods, list_of_guesses, label_batch, '', show=True)

        guesses = finetuner_guess.detach().cpu().numpy()
        labels = label_batch[:, :-1].detach().cpu().numpy()

        from scipy.stats import kstest
        import matplotlib.pyplot as plt
        import numpy as np

        score_0 = kstest(guesses[:, 0], labels[:, 0])
        score_2 = kstest(guesses[:, 2], labels[:, 2])
        score_26 = kstest(guesses[:, 26], labels[:, 26])
        score_48 = kstest(guesses[:, 48], labels[:, 48])

        print("score_0", score_0)
        print("score_2", score_2)
        print("score_26", score_26)
        print("score_48", score_48)

        scores = []
        for i in range(guesses.shape[1]):
            scores.append(kstest(guesses[:, i], labels[:, i])[0])
        scores = np.sort(np.array(scores))

        print("Model:", model_name, np.mean(scores[int(scores.shape[0]*.75):]))

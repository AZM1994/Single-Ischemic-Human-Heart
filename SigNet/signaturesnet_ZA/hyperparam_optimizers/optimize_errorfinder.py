import gc
import os
import sys

from skopt.space import Real, Integer, Categorical
import torch
import wandb

from signaturesnet.HyperParameterOptimizer import GaussianProcessSearch, ParallelSearcher
from signaturesnet.trainers import train_errorfinder
from signaturesnet.utilities.io import read_config

experiment_id = "errorfinder_lagrange"

# Training params
# batch_sizes = Integer(name='batch_size', low=20, high=500)
# learning_rates = Real(name='lr', low=1e-5, high=1e-3)
# neurons_pos = Integer(name='num_neurons_pos', low=50, high=1200)
# layers_pos = Integer(name='num_hidden_layers_pos', low=1, high=6)
# neurons_neg = Integer(name='num_neurons_neg', low=50, high=1200)
# layers_neg = Integer(name='num_hidden_layers_neg', low=1, high=6)

# Loss params
lagrange_base = Real(name="lagrange_base", low=0.01, high=1.0)
lagrange_high_error_sigs = Real(name="lagrange_high_error_sigs", low=0.01, high=1.0)
# lagrange_pnorm_param = Real(name="lagrange_pnorm", low=1e3, high=1e5)
# lagrange_smalltozero_param = Real(name="lagrange_smalltozero", low=0.1, high=10)
# pnorm_order_param = Integer(name="pnorm_order", low=3, high=9)

input_file = "search_results_" + experiment_id + ".csv"
output_file = "search_results_" + experiment_id + "_2.csv"

if __name__ == "__main__":
# Select training device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    config = read_config(path="../configs/errorfinder/errorfinder.yaml")

    # Define hyperparameters to train
    search_space = [lagrange_base, lagrange_high_error_sigs]

    # search_space = [batch_sizes, learning_rates, neurons_pos,
                    # layers_pos, neurons_neg, layers_neg,
                    # pnorm_order_param, lagrange_smalltozero_param,
                    # lagrange_pnorm_param, lagrange_base, lagrange_high_error_sigs]
    fixed_space = {"plot": True}

    def objective(**kwargs):
        loss_params = {
            "lagrange_base": kwargs["lagrange_base"],
            "lagrange_high_error_sigs": kwargs["lagrange_high_error_sigs"],
            "lagrange_pnorm": 33434.0,
            "lagrange_smalltozero":  0.276,
            "pnorm_order": 5.0,
        }

        run_config = {#"batch_size": config["batch_size"],
        #             "lr": config["lr"],
        #             "num_neurons_pos": config["num_neurons_pos"],
        #             "num_hidden_layers_pos": config["num_hidden_layers_pos"],
        #             "num_neurons_neg": config["num_neurons_neg"],
        #             "num_hidden_layers_neg": config["num_hidden_layers_neg"],
                    "loss_params": loss_params,
                    "models_dir": "../../trained_models/exp_final"
                    }
        config.update(run_config)
        print(config)
        val = train_errorfinder(config=config, data_folder="../../data")
        return val

    gp_search = GaussianProcessSearch(search_space=search_space,
                                      fixed_space=fixed_space,
                                      evaluator=objective,
                                      input_file=input_file,  # Use None to start from zero
                                      output_file=output_file)  # Store tested points
    gp_search.init_session()
    best_metaparams, best_val = gp_search.get_maximum(
        n_calls=1000,
        n_random_starts=2,
        noise=0.01,
        verbose=True,
        plot_results=True)

    print("best_metaparams:", best_metaparams)
    print("best_val:", best_val)

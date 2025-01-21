import os
import sys

from skopt.space import Real, Integer, Categorical
import torch
import wandb

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.io import read_config, read_data
from trainers import train_finetuner
from HyperParameterOptimizer.src.gaussian_process_search import GaussianProcessSearch

assert len(sys.argv[1]) == 2, "Usage: python optimize_finetuner low"
large_low = str(sys.argv[1])
assert large_low in ["low", "large"], "Argument can only be large or low"

experiment_id = "fixed_bayesian_finetuner_oversample_%s"%large_low
num_classes = 72

batch_sizes = Integer(name='batch_size', low=20, high=500)
learning_rates = Real(name='lr', low=1e-6, high=5e-3)
num_neurons = Integer(name='num_neurons', low=50, high=1000)
num_hidden_layers = Integer(name='num_hidden_layers', low=1, high=6)

input_file = None  # Use None to start from zero
output_file = "search_results_" + experiment_id + ".csv"

if __name__ == "__main__":
    # Read base config
    config = read_config(path="../configs/finetuner/finetuner_%s.yaml"%large_low)

    # Define hyperparameters to train
    search_space = [batch_sizes, learning_rates, num_hidden_layers, num_neurons]
    fixed_space = {"plot": True}

    def objective(batch_size,
                  lr,
                  num_hidden_layers,
                  num_neurons,
                  plot=False):
        run_config = {
            "model_id": experiment_id + "_run",
            "wandb_project_id": experiment_id,
            "batch_size": batch_size,
            "lr": lr,
            "num_hidden_layers": num_hidden_layers,
            "num_neurons": num_neurons,
            }
        config.update(run_config)
        val = train_finetuner(config=config, data_folder="../../data/", name=str(num_neurons))
        return val

    # Start optimization
    gp_search = GaussianProcessSearch(search_space=search_space,
                                      fixed_space=fixed_space,
                                      evaluator=objective,
                                      input_file=input_file,
                                      output_file=output_file)  # Store tested points
    gp_search.init_session()
    best_metaparams, best_val = gp_search.get_maximum(
        n_calls=60,
        n_random_starts=30,
        noise=0.01,
        verbose=True,
        plot_results=True)

    print("best_metaparams:", best_metaparams)
    print("best_val:", best_val)

import os
import sys

from skopt.space import Real, Integer, Categorical
import torch
import wandb

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.io import read_config
from trainers import train_generator, log_results
from HyperParameterOptimizer.src.gaussian_process_search import GaussianProcessSearch

experiment_id = "bayesian_generator"

lagrange_param = Real(name='lagrange_param', low=0.00001, high=0.9)
batch_size = Integer(name='batch_size', low=10, high=500)
learning_rate_encoder = Real(name='learning_rate_encoder', low=1e-5, high=5e-3)
learning_rate_decoder = Real(name='learning_rate_decoder', low=1e-5, high=5e-3)
latent_dim = Integer(name='latent_dim', low=10, high=200)
num_hidden_layers = Integer(name='num_hidden_layers', low=1, high=8)

input_file = None  # Use None to start from zero
output_file = "search_results_" + experiment_id + ".csv"

if __name__ == "__main__":
    # Read base config
    config = read_config(path="../configs/generator/generator.yaml")
    

    # Define hyperparameters to train
    search_space = [lagrange_param, batch_size, learning_rate_encoder, learning_rate_decoder, latent_dim, num_hidden_layers]
    fixed_space = {"plot": True}

    def objective(lagrange_param,
                  batch_size,
                  learning_rate_encoder,
                  learning_rate_decoder,
                  latent_dim,
                  num_hidden_layers,
                  plot=False):

        run_config = {
            "model_id": experiment_id + "_run",
            "wandb_project_id": experiment_id,
            # Parameters
            "lagrange_param": lagrange_param,
            "batch_size": batch_size,
            "learning_rate_encoder": learning_rate_encoder,
            "learning_rate_decoder": learning_rate_decoder,
            "latent_dim": latent_dim,
            "num_hidden_layers": num_hidden_layers,
            }

        config.update(run_config)
        train_DQ99R = train_generator(config=config, data_folder="../../data/")
        return -train_DQ99R

    # Start optimization
    gp_search = GaussianProcessSearch(search_space=search_space,
                                      fixed_space=fixed_space,
                                      evaluator=objective,
                                      input_file=input_file,
                                      output_file=output_file)  # Store tested points
    gp_search.init_session()
    best_metaparams, best_val = gp_search.get_maximum(
        n_calls=200,
        n_random_starts=50,
        noise=0.02,
        verbose=True,
        plot_results=True)

    print("best_metaparams:", best_metaparams)
    print("best_val:", best_val)

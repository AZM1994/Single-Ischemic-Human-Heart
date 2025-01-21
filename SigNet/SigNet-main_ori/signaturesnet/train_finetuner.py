import os
import sys

from argparse import ArgumentParser
import torch
import wandb

from signaturesnet.trainers import train_finetuner
from signaturesnet.utilities.io import update_dict, read_config, write_result

# DEFAULT_CONFIG_FILE = ["configs/finetuner/finetuner_low.yaml"]
DEFAULT_CONFIG_FILE = ["configs/finetuner/finetuner_large.yaml"]

if __name__ == "__main__":
    # Parse command-line arguments
    parser = ArgumentParser()
    parser.add_argument(
        '--config_file', action='store', nargs=1, type=str, required=False, default=DEFAULT_CONFIG_FILE,
        help=f'Path to yaml file containing all needed parameters.\nThey can be overwritten by command-line arguments'
    )

    parser.add_argument(
        '--model_id', action='store', nargs=1, type=str, required=False,
        help=f'Unique id given to the trained model.'
    )

    # Data args
    parser.add_argument(
        '--source', action='store', nargs=1, type=str, required=False,
        help='Data source.'
    )
    parser.add_argument(
        '--network_type', action='store', nargs=1, type=str, required=False,
        help='Network type: either random or realistic.'
    )
    parser.add_argument(
        '--sigmoid_params', action='store', nargs=1, type=list, required=False,
        help='Sigmoid parameters for normalization.'
    )
    parser.add_argument(
        '--load_data', action='store', nargs=1, type=str, required=False,
        help='Whether to load or generate data for the training'
    )

    # Train args
    parser.add_argument(
        '--batch_size', action='store', nargs=1, type=int, required=False,
        help='Training batch size.'
    )
    parser.add_argument(
        '--lr', action='store', nargs=1, type=float, required=False,
        help='Learning rate.'
    )
    
    # Model args
    parser.add_argument(
        '--num_neurons', action='store', nargs=1, type=int, required=False,
        help='Neurons of the positive error network.'
    )
    parser.add_argument(
        '--num_hidden_layers', action='store', nargs=1, type=int, required=False,
        help='Layers of the positive error network.'
    )

    _args = parser.parse_args()

    # Read & config
    config = read_config(path=getattr(_args, "config_file")[0])
    config = update_dict(config=config, args=_args)
    
    print("Using config file:", getattr(_args, "config_file")[0])
    print("Using config:", config)
    score = train_finetuner(config=config)
    write_result(score, "../tmp/finetuner_%s_score_%s.txt"%(config["network_type"],config["model_id"]))
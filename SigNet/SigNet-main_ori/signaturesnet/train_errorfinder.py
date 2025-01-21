from argparse import ArgumentParser
import os
import sys

import torch
import wandb

from signaturesnet.trainers import train_errorfinder
from signaturesnet.utilities.io import update_dict, read_config, write_result

DEFAULT_CONFIG_FILE = ["configs/errorfinder/errorfinder.yaml"]
torch.manual_seed(123)

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
    parser.add_argument(
        '--source', action='store', nargs=1, type=str, required=False,
        help=f'Data source to use in the training.'
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
    parser.add_argument(
        '--sigmoid_params', action='store', nargs=1, type=list, required=False,
        help='Sigmoid parameters for normalization.'
    )
    # Model args
    parser.add_argument(
        '--num_neurons_pos', action='store', nargs=1, type=int, required=False,
        help='Neurons of the positive error network.'
    )
    parser.add_argument(
        '--num_hidden_layers_pos', action='store', nargs=1, type=int, required=False,
        help='Layers of the positive error network.'
    )
    parser.add_argument(
        '--num_neurons_neg', action='store', nargs=1, type=int, required=False,
        help='Neurons of the negative error network.'
    )
    parser.add_argument(
        '--num_hidden_layers_neg', action='store', nargs=1, type=int, required=False,
        help='Layers of the negative error network.'
    )

    # Loss args
    parser.add_argument(
        '--lagrange_base', action='store', nargs=1, type=float, required=False,
        help='Weight of misclassification per base signatures.'
    )
    parser.add_argument(
        '--lagrange_high_error_sigs', action='store', nargs=1, type=float, required=False,
        help='Weight of misclassification per high error signatures.'
    )
    parser.add_argument(
        '--lagrange_pnorm', action='store', nargs=1, type=float, required=False,
        help='Weight of big errors in the loss.'
    )
    parser.add_argument(
        '--lagrange_smalltozero', action='store', nargs=1, type=float, required=False,
        help='Weight of not sending small values to zero in the loss.'
    )
    parser.add_argument(
        '--pnorm_order', action='store', nargs=1, type=int, required=False,
        help='Norm used in the loss.'
    )
    _args = parser.parse_args()

    # Read & config
    config = read_config(path=getattr(_args, "config_file")[0])
    config = update_dict(config=config, args=_args)
    config["loss_params"] = update_dict(config=config["loss_params"], args=_args)
    
    print("Using config:", config)
    score = train_errorfinder(config=config)
    # write_result(score, "../tmp/errorfinder_score_%s_%s.txt"%(config["source"],config["model_id"]))
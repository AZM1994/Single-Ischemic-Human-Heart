import os
import logging

#import numpy as np
from argparse import ArgumentParser
import pandas as pd
import torch

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.modules.signet_module import SigNet
from signaturesnet.models import Generator
from signaturesnet.models import Classifier
from signaturesnet.utilities.io import read_model, tensor_to_csv
from signaturesnet.utilities.normalize_data import normalize_data

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        'task', help=f'Must be one of: [refitter, generator, detector]. Depending on the solution you are interested in.'
    )

    parser.add_argument(
        '--input_data', action='store', nargs=1, type=str, required=False,
        help=f'Path to the input data to be analyzed. By default it will use PCAWG dataset'
    )

    parser.add_argument(
        '--normalization', action='store', nargs=1, type=str, required=False, default=[None],
        help=f'The kind of normalization to be applied to the data. Should be either "None" (default), "exome", "genome" or a path to a file with the oppportunities.'
    )

    parser.add_argument(
        '--output_path', action='store', nargs=1, type=str, required=False, default="signet_output",
        help=f"Name of this inference's results"
    )

    parser.add_argument(
        '--plot_figs', action='store', nargs=1, type=bool, required=False, default=[False],
        help=f'Boolean. Whether to compute plots for the output. Default: "False".'
    )

    parser.add_argument(
        '--n_points', action='store', nargs=1, type=str, required=False, default=[1000],
        help=f'[ONLY FOR GENERATOR] Number of points to be generated.'
    )

    
    args = parser.parse_args()
    return args


def run_refitter(args):
    # Read data
    mutations = pd.read_csv(args.input_data[0], header=0, index_col=0)

    # Load & Run signet
    signet = SigNet(opportunities_name_or_path=args.normalization[0])
    results = signet(mutation_dataset=mutations)

    # Store results
    results.save(path=args.output_path[0])

    # Plot figures
    if args.plot_figs:
        results.plot_results(save=True)

def run_generator(args):
    # Read model
    logging.info("Loading model...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    generator = read_model(os.path.join(TRAINED_MODELS, "generator"), device)
    # generator = Generator().to(device)

    # Generate data
    logging.info("Generating data...")
    labels = generator.generate(batch_size=int(args.n_points[0]))
    
    # Store results
    logging.info(f"Saving generated data to {args.output_path[0]}")
    filepath = os.path.join(args.output_path[0], "generated_data.csv")
    tensor_to_csv(labels, filepath)
    
def run_detector(args):
    # Read model
    logging.info("Loading model...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    detector = read_model(os.path.join(TRAINED_MODELS, "detector"), device)

    logging.info("Reading data...")
    mutation_dataset = pd.read_csv(args.input_data[0], header=0, index_col=0)
    
    logging.info("Preprocessing input data...")
    mutation_order = pd.read_excel(os.path.join(DATA, "mutation_type_order.xlsx"))
    mutation_dataset = mutation_dataset[list(mutation_order['Type'])]
    sample_names = mutation_dataset.index  # NOTE(Oleguer): Should we be using this in the output to notify the order?
    mutation_vec = torch.tensor(mutation_dataset.values, dtype=torch.float, device=device)

    logging.info("Normalizing input data...")
    if args.normalization[0] is not None:
        print(args.normalization)
        mutation_vec = normalize_data(mutation_vec, args.normalization[0])
    sums = torch.sum(mutation_vec, dim=1).reshape(-1, 1)
    normalized_mutation_vec = mutation_vec / sums

    logging.info("Detecting out-of-train-distribution points...")
    num_mutations = torch.sum(mutation_vec, dim=1)
    classification = detector(mutation_dist=normalized_mutation_vec,
                              num_mut=num_mutations.reshape(-1, 1)).view(-1)

    # Store results
    logging.info(f"Saving classification results into {args.output_path[0]}")
    filepath = os.path.join(args.output_path[0], "classification_results.csv")  #TODO(Oleguer): Should we add an index or something?
    tensor_to_csv(classification, filepath)

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    VALID_TASKS = ['refitter', 'generator', 'detector']
    assert args.task in VALID_TASKS, f"Task must be one of {VALID_TASKS}. You provided {args.task}"
    print(args.input_data)

    if args.task == "refitter":
        assert args.input_data[0] is not None, "You must provide an input data to run the refitter. --input_data=<your_path>"
        run_refitter(args)
    elif args.task == "generator":
        assert args.n_points[0] is not None, "You must provide the number of points to generate. --n_points=1000"
        run_generator(args)
    elif args.task == "detector":
        assert args.input_data[0] is not None, "You must provide an input data to run the refitter. --input_data=<your_path>"
        run_detector(args)


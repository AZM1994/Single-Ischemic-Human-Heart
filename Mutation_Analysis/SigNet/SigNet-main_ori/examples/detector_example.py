import os
import logging

from argparse import ArgumentParser
import torch

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.utilities.io import csv_to_tensor, read_model


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        '--input_data', action='store', nargs=1, type=str, required=False, default=[DATA + "/datasets/example_input.csv"],
        help=f'Path to the input data to be analyzed. By default it will use PCAWG dataset'
    )
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()
    
    # Read data
    mutation_vectors = csv_to_tensor(args.input_data[0], header=0, index_col=0)
    
    # Get the number of mutations for each mutational vector and normalize the vector
    num_muts = torch.sum(mutation_vectors, dim=1).reshape(-1, 1)
    normalized_mutation_vectors = mutation_vectors / num_muts

    # Run classification
    detector = read_model(os.path.join(TRAINED_MODELS, "detector"))
    classification_guess = detector(mutation_dist=normalized_mutation_vectors,
                                    num_mut=num_muts)

    # Classify results
    classification_cutoff = 0.5
    classification_results = (classification_guess >= classification_cutoff).to(torch.int64)
    
    print("Classifications:\n", classification_results)
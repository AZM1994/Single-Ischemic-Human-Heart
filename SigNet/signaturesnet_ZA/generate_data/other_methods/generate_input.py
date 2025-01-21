import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from models import create_baseline_dataset
from utilities.io import read_signatures, sort_signatures

signatures_file = "../../data/data_v2.xlsx"
input_file = "../../data/real_data/PCAWG_norm.csv"
output_signatures = "PCAWG_output/signatures.csv"
output_file = "PCAWG_output/PCAWG_baseline.csv"

# Sort and save signatures file:
signatures = sort_signatures(signatures_file, output_signatures)

# Apply baseline with corresponding signatures to PCAWG data:
create_baseline_dataset(input_file, output_file, signatures_file, which_baseline="slow")
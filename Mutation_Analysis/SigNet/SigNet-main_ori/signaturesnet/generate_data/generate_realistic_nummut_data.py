import itertools

import torch

from signaturesnet.utilities.io import tensor_to_csv
from signaturesnet.utilities.temporal_io import generate_realistic_nummut_data
from signaturesnet import DATA

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

if __name__ == "__main__":
    splits = ["train", "val", "test"]
    large_or_low = ["large", "low"]
    n_samples = {"train": 10_000, "val": 1_000, "test": 1_000}

    for split, large_low in itertools.product(splits, large_or_low):
        data = generate_realistic_nummut_data(n_samples=n_samples[split],
                                              split=split,
                                              large_or_low=large_low,
                                              device=DEVICE,)
        tensor_to_csv(data.inputs, DATA + f"/realistic_nummuts_data/{split}_{large_low}_input.csv")
        tensor_to_csv(data.labels, DATA + f"/realistic_nummuts_data/{split}_{large_low}_label.csv")
        tensor_to_csv(data.prev_guess, DATA + f"/realistic_nummuts_data/{split}_{large_low}_baseline.csv")


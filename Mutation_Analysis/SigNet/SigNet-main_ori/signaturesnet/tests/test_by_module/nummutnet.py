import os

import torch

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.utilities.io import read_model

if __name__=="__main__":
    model = read_model(os.path.join(TRAINED_MODELS, "nummutnet"))
    print(model)
    nummuts = model.get_nummuts(torch.rand(5, 72))
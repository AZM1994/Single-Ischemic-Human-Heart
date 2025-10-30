import os

from signaturesnet import DATA, TRAINED_MODELS
from signaturesnet.utilities.io import read_model, read_signatures

if __name__=="__main__":
    signatures = read_signatures(os.path.join(DATA, "data.xlsx"))
    generator = read_model(os.path.join(TRAINED_MODELS, "generator"))
    nummutnet = read_model(os.path.join(TRAINED_MODELS, "nummutnet"))

    mutvec, weights = generator.generate_weights_and_samples(n_samples=10,
                                                             signatures=signatures, 
                                                             realistic_number_of_mutations=True,
                                                             nummutnet=nummutnet,)

    print("mutvec\n", mutvec)
    print("weights\n", weights)
import torch
import pandas as pd

from signaturesnet.utilities.io import csv_to_tensor, read_signatures, tensor_to_csv


# Function to sample:
def sample_from_sig(signature, num_mut, normalize):
    sample = signature
    c = torch.distributions.categorical.Categorical(probs=signature)
    sample_shape = torch.Size([num_mut, ])
    samples = c.sample(sample_shape=sample_shape).type(torch.float32)
    n_bins = signature.shape[0]
    sample = torch.histc(samples, bins=n_bins, min=0, max=n_bins - 1)
    if normalize:
        sample = sample/float(num_mut)
    return sample

# Read signature weights data:
real_weights = csv_to_tensor("../data/real_data/sigprofiler_not_norm_PCAWG.csv", device='cpu', header=0, index_col=0)
real_weights = torch.cat([real_weights, torch.zeros(real_weights.size(0), 7).to(real_weights)], dim=1) 

# Read signatures
signatures = read_signatures(
    "../data/data.xlsx", mutation_type_order="../data/mutation_type_order.xlsx")

# List of numbers of mutations
num_muts = list(range(1,25))

# Create input-labels pairs
input_batch = torch.empty((real_weights.shape[0]*len(num_muts), 96))
labels_batch = torch.empty((real_weights.shape[0]*len(num_muts), 73))
for i, num_mut in enumerate(num_muts):
    for j in range(len(real_weights)):
        # Compute resulting signature
        signature = torch.einsum("ij,j->i", (signatures, real_weights[j]))

        # Sample the mutational vector from the given distribution
        sample = sample_from_sig(signature=signature,
                                num_mut=int(num_mut),
                                normalize=True)
        # Store Values
        input_batch[i*len(real_weights)+j, :] = sample
        labels_batch[i*len(real_weights)+j, :] = torch.cat([real_weights[j, :], torch.tensor([float(num_mut)])])

# Convert to dataframe and add colnames and rownames:
mutation_order = pd.read_excel("../data/mutation_type_order.xlsx")
input_batch = pd.DataFrame(input_batch.detach().numpy())
input_batch.columns = mutation_order['Type']
input_batch.index = ['sample_'+str(i) for i in range(len(input_batch))]

# Save results
input_batch.to_csv('../../data/exp_superlow_nummut/inputs.csv')
tensor_to_csv(labels_batch, '../../data/exp_superlow_nummut/labels.csv')

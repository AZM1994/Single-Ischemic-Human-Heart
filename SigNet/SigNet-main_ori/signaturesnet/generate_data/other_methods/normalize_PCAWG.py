import os
import sys
import pandas as pd
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.io import read_signatures
from utilities.normalize_data import create_opportunities


data = pd.read_csv("../../data/real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv", index_col=0)
cols = data.columns.tolist()
rows = data.index.tolist()

num_muts = data.sum(axis=1, numeric_only=True)
data = torch.tensor(data.values, dtype=torch.float)

signatures = read_signatures("../../data/sigprofiler_sigs.xlsx")
abundances = torch.Tensor(create_opportunities('../../data/data_donors/3mer_WG_hg37.txt'))

w = torch.zeros(data.size())
for i in range(data.shape[1]):
    for j in range(data.shape[0]):
        w[j,i] = torch.sum(torch.div(data[j,i]*signatures[:,i], abundances))

w = torch.div(w,torch.sum(w,dim=1).reshape(-1,1))
w = torch.multiply(w,torch.tensor(num_muts).reshape(-1,1))

df = pd.DataFrame(w.detach().numpy())
df.columns = cols
row_names =rows
df.index = row_names
df.to_csv("../../data/real_data/sigprofiler_normalized_PCAWG.csv", header=True, index=True)



#### Normalize to sum 1 without taking into consideration abundances:

w = torch.div(data,torch.sum(data,dim=1).reshape(-1,1))
df = pd.DataFrame(w.detach().numpy())
df.columns = cols
row_names =rows
df.index = row_names
df.to_csv("../../data/real_data/sigprofiler_not_norm_PCAWG.csv", header=True, index=True)

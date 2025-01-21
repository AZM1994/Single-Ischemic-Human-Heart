import sys
import os

import copy
import math
import numpy as np
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from signaturesnet import DATA
from signaturesnet.utilities.plotting import plot_correlation_matrix
from signaturesnet.utilities.io import sort_signatures

class OverSampler:
    def __init__(self, data):
        self.data = data
        self.presence = self._get_presence(data)

    def _get_presence(self, data, offset = 0.05):
        presence = torch.sum(data, dim=0)
        return presence / torch.sum(presence)

    def get_frequencies(self, points):
        frequencies = torch.einsum("i,bi->b", self.presence, points)
        return frequencies

    def get_oversampled_set(self, rarity=0.1, n_repetitions=5):
        """Take the rarity% and repeats it n_repetitions
        """
        oversampled = copy.deepcopy(self.data)
        frequencies = self.get_frequencies(oversampled)
        sorted_frequencies, indices = torch.sort(frequencies)
        to_repeat = oversampled[indices[:int(rarity*indices.size(0))], ...]
        to_append = [to_repeat]*n_repetitions
        oversampled = torch.cat([oversampled] + to_append, dim=0)
        return oversampled[np.random.permutation(oversampled.size(0)), ...]


class CancerTypeOverSampler:
    def __init__(self, data, cancer_types):
        self.data = data
        self.cancer_types = cancer_types
        self.counts = torch.bincount(cancer_types.to(torch.int))
        # print(cancer_types)
        # print(self.counts)

    def get_oversampled_set(self, rarity=0.1, n_repetitions=5):
        """Take the rarity% and repeats it n_repetitions
        """
        oversampled = copy.deepcopy(self.data)
        uniform_count = int(self.cancer_types.size(0)/self.counts.size(0))
        for i in range(self.cancer_types.size(0)):
            count = self.counts[int(self.cancer_types[i].item())]
            times_to_repeat = min(int(torch.round(uniform_count/count)), 10)  # dont repeat more than 10 times
            # print(i, int(self.cancer_types[i].item()), count, times_to_repeat)
            to_append = [self.data[i, :].unsqueeze(0)]*times_to_repeat
            oversampled = torch.cat([oversampled] + to_append, dim=0)
        return oversampled[np.random.permutation(oversampled.size(0)), ...]

    def get_even_set(self):
        """Oversample to have the same number of samples for each cancer type
        """
        oversampled = torch.tensor([])
        max_count = torch.max(self.counts)
        for i in range(self.cancer_types.size(0)):
            count = self.counts[int(self.cancer_types[i].item())]
            times_to_repeat = int(max_count/count)
            # print(i, int(self.cancer_types[i].item()), count, times_to_repeat)
            to_append = [self.data[i, :].unsqueeze(0)]*times_to_repeat
            oversampled = torch.cat([oversampled] + to_append, dim=0)
        return oversampled[np.random.permutation(oversampled.size(0)), ...]

    def get_N_oversampled_set(self, N_samples):
        """Oversample to have at least N number of samples for each cancer type
        """
        oversampled = torch.tensor([])
        for i in range(self.cancer_types.size(0)):
            count = self.counts[int(self.cancer_types[i].item())]
            times_to_repeat = math.ceil(N_samples/int(count))
            # print(i, int(self.cancer_types[i].item()), count, times_to_repeat)
            to_append = [self.data[i, :].unsqueeze(0)]*times_to_repeat
            oversampled = torch.cat([oversampled] + to_append, dim=0)
        return oversampled[np.random.permutation(oversampled.size(0)), ...]


if __name__=="__main__":
    import os
    import sys

    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from utilities.io import csv_to_tensor, csv_to_pandas

    real_data = csv_to_pandas(DATA + "real_data/sigprofiler_not_norm_PCAWG.csv",
                              device="cpu", header=0, index_col=0,
                              type_df=DATA + "real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv")
    # print(real_data)
    # real_data_g = real_data.groupby('cancer_type').sample(frac=1)
    # print(real_data_g.groupby('cancer_type').head(3))
    # print(real_data_g.groupby('cancer_type').tail(3))
    # print(real_data_g['cancer_type'][-1])
    # print(real_data.shape[0])

    print(real_data)

    real_data = torch.tensor(real_data.values, dtype=torch.float)
    real_data, cancer_types = real_data[:, :-1], real_data[:, -1]

    real_data = real_data/torch.sum(real_data, axis=1).reshape(-1, 1)
    real_data = torch.cat([real_data, torch.zeros(real_data.size(0), 7).to(real_data)], dim=1)

    oversampler = CancerTypeOverSampler(real_data, cancer_types)
    new_set = oversampler.get_N_oversampled_set(500)

    print(real_data.size())
    print(new_set.size())

    print(sum(real_data[:,2]>0))
    print(sum(new_set[:,2]>0))
    print(sum(real_data[:,0]>0))
    print(sum(new_set[:,0]>0))

    data_folder = DATA
    signatures = sort_signatures(file=data_folder + "data.xlsx",
                                 mutation_type_order=data_folder + "mutation_type_order.xlsx")
    corrMatrix = plot_correlation_matrix(data=new_set, signatures=signatures)
    corrMatrix = plot_correlation_matrix(data=real_data, signatures=signatures)
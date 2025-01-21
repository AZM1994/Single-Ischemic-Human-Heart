import os
import sys
import logging

import numpy as np
import random
import torch
from tqdm import tqdm

from signaturesnet.utilities.io import csv_to_tensor, read_model

class DataGenerator:

    def __init__(self,
                 signatures,
                 seed=None,
                 normalize=True,
                 shuffle=True):
        self.signatures = signatures
        self.total_signatures = signatures.shape[1]
        self.normalize = normalize
        self.shuffle = shuffle
        self._set_all_seeds(seed)

        logging.info("Total signatures: %s"%self.total_signatures)


    def _set_all_seeds(self, seed):
        if seed is None:
            return
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.deterministic = True


    def _get_nummuts(self, split, large_or_low, size):
        assert split in ["train", "valid", "test"]
        assert large_or_low in ["large", "low", "superlow"]

        # For train and validation we sample from ranges
        if split == "train" or split == "valid":
            if large_or_low == "superlow": 
                range_muts = [1, 5, 10, 15, 20, 25]
                n_ranges = len(range_muts) - 1
                num_samples = int(size/n_ranges)
                ind_range_muts = [0]*num_samples + [1]*num_samples + [2]*num_samples + \
                    [3]*num_samples + [4]*(num_samples + size%n_ranges)

            elif large_or_low == 'low':
                range_muts = [15, 50, 100, 250, 500, 1000, 5000, 1e4, 1e5]
                n_ranges = len(range_muts) - 1
                num_samples = int(size/n_ranges)
                ind_range_muts = [0]*num_samples + [1]*num_samples + [2]*num_samples + \
                    [3]*num_samples + [4]*num_samples + [5]*num_samples + [6]*num_samples +\
                    [7]*(num_samples + size%n_ranges)
            
            elif large_or_low == 'large':
                range_muts = [1e3, 5e3, 1e4, 5e4, 1e5, 5e5]
                n_ranges = len(range_muts) - 1
                num_samples = int(size/n_ranges)
                # ind_range_muts = [0]*num_samples + [1]*num_samples + \
                #     [2]*num_samples + [3]*num_samples + [-1]*(num_samples + size%8)
                ind_range_muts = [0]*num_samples + [1]*num_samples + \
                    [2]*num_samples + [3]*num_samples + [-1]*(num_samples)

                
#                print("n_ranges, num_samples, ind_range_muts", n_ranges, num_samples, ind_range_muts, len(ind_range_muts), len([0]*num_samples))
                # -1 means real distribution

            nummuts = [np.random.randint(range_muts[ind_range_muts[i]], range_muts[ind_range_muts[i] + 1])\
                        if ind_range_muts[i] != -1 else -1 for i in range(len(ind_range_muts))]
            return nummuts

        # For testing we use exact number of mutations
        num_samples = int(size/10.0)
        num_muts = [25]*num_samples + [50]*num_samples + [100]*num_samples +\
                   [250]*num_samples + [500]*num_samples + [1e3]*num_samples +\
                   [5e3]*num_samples + [1e4]*num_samples + [5e4]*num_samples +\
                   [1e5]*(num_samples + size%10)
        return num_muts


    def _sample_from_sig(self, signature, num_mut):
        sample = signature
        if num_mut > 0 and num_mut < 1e5:  # only do sampling for num_mut < 1e5
            c = torch.distributions.categorical.Categorical(probs=signature)
            sample_shape = torch.Size([num_mut, ])
            samples = c.sample(sample_shape=sample_shape).type(torch.float32)
            n_bins = signature.shape[0]
            sample = torch.histc(samples, bins=n_bins, min=0, max=n_bins - 1)
            if self.normalize:
                sample = sample/float(num_mut)
        else:
            sample = sample/torch.sum(sample)
        return sample


    def make_input(self, labels, split, large_low, augmentation=1, nummuts=None):
        """Create a labelled dataset of mutation vectors
           from a tensor of labels.
           Returns:
            inputs (mutational vector)
            labels (including an appended column with the number of mutations used)
        """
        labels = torch.cat([labels]*augmentation, dim = 0).to(torch.float32)
        if nummuts is None:
            nummuts = self._get_nummuts(split, large_low, size=labels.shape[0])
        else:
            nummut_means = torch.cat([nummuts]*augmentation, dim=0)
            nummut_stds = nummut_means/10.0
            nummuts = torch.clip(torch.distributions.normal.Normal(nummut_means, nummut_stds).sample().type(torch.int64), 10, 1e5)

        input_batch = torch.empty((labels.shape[0], 96))
        labels_batch = torch.empty((labels.shape[0], self.total_signatures + 1))

        logging.info("Sampling from labels...")
        for i, num_mut in tqdm(enumerate(nummuts)):
            # Compute resulting signature
            signature = torch.einsum("ij,j->i", (self.signatures, labels[i]))

            # Sample the mutational vector from the given distribution
            sample = self._sample_from_sig(signature=signature,
                                           num_mut=int(num_mut))
            
            # Store Values
            # For the real distribution we say we have more than 1e5 mutations
            num_mut = float(num_mut) if num_mut != -1 else float(np.random.randint(1e5, 1e6))
            input_batch[i, :] = sample
            labels_batch[i, :] = torch.cat([labels[i, :], torch.tensor([float(num_mut)])])
        return input_batch, labels_batch


    def make_random_set(self,
                        split,
                        large_low,
                        num_samples,
                        min_n_signatures=1,
                        max_n_signatures=10,
                        min_weight=0.1):
        """Create a labelled dataset of mutation vectors
           randomly combining the signatures.
        """
        assert(min_n_signatures > 0 and max_n_signatures <= self.total_signatures)

        labels = torch.empty((num_samples, self.total_signatures))

        logging.info("Generating random labels...")
        for i in tqdm(range(num_samples)):
            if i == 1:
                print("here")
            # Pick the number of signatures involved
            n_signatures = np.random.randint(min_n_signatures, max_n_signatures + 1)

            # Select n_signatures
            signature_ids = torch.randperm(self.total_signatures)[:n_signatures]

            # Assign weights randomly
            weights = torch.rand(size=(n_signatures,))
            weights = weights/torch.sum(weights)
            weights[weights < min_weight] = 0
            weights = weights/torch.sum(weights)
            label = torch.zeros(self.total_signatures).scatter_(
                dim=0, index=signature_ids, src=weights)

            # Store
            labels[i, :] = label
        return self.make_input(labels, split, large_low)


    def make_realistic_set(self,
                           generator_model_path,
                           split,
                           large_low,
                           num_samples,
                           std = 1.0):
        """Create a labelled dataset of mutation vectors
           from the generator output weights.
        """
        generator = read_model(generator_model_path, device="cpu")
        logging.info("Generating realistic labels...")
        labels = generator.generate(num_samples, std = std)
        return self.make_input(labels, split, large_low)
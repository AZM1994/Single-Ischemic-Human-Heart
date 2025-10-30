import torch
from torch.utils.data import Dataset

class GeneratorData(Dataset):
    def __init__(self, inputs=None, cancer_types=None):
        self.inputs = inputs
        self.cancer_types = cancer_types

    def __len__(self):
        return self.inputs.shape[0]

    def __getitem__(self, i):
        return self.inputs[i]
    
    def to(self, device):
        self.inputs = self.inputs.to(device)
        self.cancer_types = self.cancer_types.to(device)
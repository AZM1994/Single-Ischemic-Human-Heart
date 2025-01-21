import torch
from torch.utils.data import Dataset


class DataPartitions(Dataset):
    def __init__(self, inputs=None, labels=None, prev_guess=None, num_mut=None, extract_nummut=True):
        """Encapsulates input, label and previous step guess of a dataset.
        
        Args:
            inputs (Torch.Tensor, optional): Input to current ANN. Defaults to None.
            labels (torch.Tensor, optional): Labels of current inputs. Defaults to None.
            prev_guess (Guess of the previous model, optional): If training finetuner, 
                this is the baseline output, if training the errorfinder this is the
                finetuner output. Defaults to None.
        """
        self.inputs = inputs
        self.labels = labels
        self.prev_guess = prev_guess
        self.classification = None
        self.num_mut = num_mut if num_mut is not None else None
        if extract_nummut and labels is not None and num_mut is None:
            self.num_mut = labels[:, -1].reshape(-1,1)
            self.labels = labels[:, :-1]

    def append(self, other):
        def cat(a, b):
            if a is not None and b is not None:
                return torch.cat([a, b], dim=0) 
        self.inputs = cat(self.inputs, other.inputs)
        self.labels = cat(self.labels, other.labels)
        self.prev_guess = cat(self.prev_guess[:,:60], other.prev_guess[:,:60])
        self.classification = cat(self.classification, other.classification)
        self.num_mut = cat(self.num_mut, other.num_mut)

    def perm(self):
        indexes = torch.randperm(self.inputs.shape[0])
        def do_perm(a):
            if a is not None:
                return a[indexes, ...]
        self.inputs = do_perm(self.inputs)
        self.labels = do_perm(self.labels)
        self.prev_guess = do_perm(self.prev_guess)
        self.classification = do_perm(self.classification)
        self.num_mut = do_perm(self.num_mut)

    def to(self, device):
        device = torch.device(device)
        self.inputs = self.inputs.to(device)
        self.classification = self.classification.to(device) if self.classification is not None else self.labels
        self.labels = self.labels.to(device) if self.labels is not None else self.labels
        self.prev_guess = self.prev_guess.to(device) if self.prev_guess is not None else self.prev_guess
        self.num_mut = self.num_mut.to(device) if self.num_mut is not None else self.num_mut

    def __len__(self):
        return self.inputs.shape[0]

    def __getitem__(self, i):
        labels_i = torch.empty(0) if self.labels is None else self.labels[i]
        prev_guess_i = torch.empty(0) if self.prev_guess is None else self.prev_guess[i]
        num_mut_i = torch.empty(0) if self.num_mut is None else self.num_mut[i]
        classification_i = torch.empty(0) if self.classification is None else self.classification[i]
        return self.inputs[i], labels_i, prev_guess_i, num_mut_i, classification_i

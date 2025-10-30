import torch
import torch.nn as nn
import torch.nn.functional as F

class NumMutNet(nn.Module):
    def __init__(self,
                 hidden_dim,
                 n_layers = 3,
                 input_dim = 72,
                 output_dim = 9):
        self.init_args = locals()
        self.init_args.pop("self")
        self.init_args.pop("__class__")
        self.init_args["model_type"] = "NumMutNet"

        super(NumMutNet, self).__init__()

        self.input = nn.Linear(input_dim, hidden_dim)
        self.dropout = nn.Dropout(0.25)
        layers = []
        for i in range(n_layers):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
        self.layers = nn.ModuleList(layers)
        self.output = nn.Linear(hidden_dim, output_dim)
        self.activation = nn.LeakyReLU(0.1)

    def forward(self, x):
        x = self.dropout(self.activation(self.input(x)))
        for layer in self.layers:
            x = self.dropout(self.activation(layer(x)))
        x = self.output(x)
        return x

    def get_nummuts(self, x):
        logits = self.forward(x)
        classes = torch.argmax(logits, dim=1)
        rand_ = torch.rand(classes.shape)
        classes_copy = classes.clone()
        classes = classes.to(torch.float)

        classes[classes_copy == 0] = rand_[classes_copy == 0] * (1.5 - 1) + 1
        classes[classes_copy == 1] = rand_[classes_copy == 1] * (2 - 1.5) + 1.5
        classes[classes_copy == 2] = rand_[classes_copy == 2] * (2.5 - 2) + 2
        classes[classes_copy == 3] = rand_[classes_copy == 3] * (3 - 2.5) + 2.5
        classes[classes_copy == 4] = rand_[classes_copy == 4] * (3.5 - 3) + 3
        classes[classes_copy == 5] = rand_[classes_copy == 5] * (4 - 3.5) + 3.5
        classes[classes_copy == 6] = rand_[classes_copy == 6] * (4.5 - 4) + 4
        classes[classes_copy == 7] = rand_[classes_copy == 7] * (5 - 4.5) + 4.5

        muts = torch.pow(10, classes).to(torch.float)
        return muts


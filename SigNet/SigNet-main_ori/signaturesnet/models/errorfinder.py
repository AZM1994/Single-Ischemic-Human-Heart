import torch
import torch.nn as nn


class SignedErrorFinder(nn.Module):
    def __init__(self,
                 num_classes=72,
                 num_hidden_layers=2,
                 num_units=400,
                 sigmoid_params = [5000, 2000]):
        super(SignedErrorFinder, self).__init__()
        self.sigmoid_params = sigmoid_params
        # Num units of the weights path
        num_units_branch_1 = int(num_units*0.9)
        # Num units of the mutations path
        num_units_branch_2 = int((num_units - num_units_branch_1)/2)
        num_units_branch_3 = num_units - num_units_branch_1 - num_units_branch_2

        # Weights path
        self.layer1_1 = nn.Linear(num_classes, num_units_branch_1)
        nn.init.kaiming_uniform_(self.layer1_1.weight)
        self.layer2_1 = nn.Linear(num_units_branch_1, num_units_branch_1)
        nn.init.kaiming_uniform_(self.layer2_1.weight)

        # Mutations path
        self.layer1_2 = nn.Linear(1, num_units_branch_2)
        nn.init.kaiming_uniform_(self.layer1_2.weight)

        # Mutations path
        self.layer1_3 = nn.Linear(1, num_units_branch_3)
        nn.init.kaiming_uniform_(self.layer1_3.weight)

        # Combined layers
        self.hidden_layers = nn.ModuleList(
            modules=[nn.Linear(num_units, num_units)
                     for _ in range(num_hidden_layers)])
        self.output_layer = nn.Linear(num_units, num_classes)
        nn.init.kaiming_uniform_(self.output_layer.weight)
        for hidden_layer in self.hidden_layers:
            nn.init.kaiming_uniform_(hidden_layer.weight)

        # Activation function
        self.activation = nn.LeakyReLU(0.1)

    def __clamp(self, x, slope=1e-2):
        return nn.LeakyReLU(slope)(1 - nn.LeakyReLU(slope)(1 - x))

    def forward(self, weights, num_mutations, classification):
        # Baseline head
        weights = self.activation(self.layer1_1(weights))
        weights = self.activation(self.layer2_1(weights))

        # Mutations head
        num_mutations = torch.log10(num_mutations)/6
        # num_mutations = nn.Sigmoid()((num_mutations-self.sigmoid_params[0])/self.sigmoid_params[1])
        num_mutations = self.activation(self.layer1_2(num_mutations))

        # Baseline head
        classification = self.activation(self.layer1_3(classification))

        # Concatenate
        comb = torch.cat([weights, num_mutations, classification], dim=1)

        # Apply shared layers
        for layer in self.hidden_layers:
            comb = self.activation(layer(comb))

        # Apply output layer
        # res = self.__clamp(self.output_layer(comb))
        res = self.output_layer(comb)
        return res


class ErrorFinder(nn.Module):

    def __init__(self,
                 num_classes=72,
                 num_hidden_layers_pos=2,
                 num_hidden_layers_neg=2,
                 num_units_pos=400,
                 num_units_neg=400):
        self.init_args = locals()
        self.init_args.pop("self")
        self.init_args.pop("__class__")
        self.init_args["model_type"] = "ErrorFinder"
        
        super(ErrorFinder, self).__init__()

        self.positive_path = SignedErrorFinder(num_classes=num_classes,
                                               num_hidden_layers=num_hidden_layers_pos,
                                               num_units=num_units_pos)

        self.negative_path = SignedErrorFinder(num_classes=num_classes,
                                               num_hidden_layers=num_hidden_layers_neg,
                                               num_units=num_units_neg)

    def forward(self, weights, num_mutations, classification):
        pred_upper = self.positive_path(weights, num_mutations, classification)
        pred_lower = self.negative_path(weights, num_mutations, classification)

        # If in eval mode, correct interval issues
        if not self.training:
            pred_upper[pred_upper < weights] = weights[pred_upper < weights]
            pred_lower[pred_lower > weights] = weights[pred_lower > weights]
            pred_lower[pred_lower < 0] = 0
            pred_upper[pred_upper > 1] = 1

        return pred_upper, pred_lower

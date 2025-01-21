import os
import sys

import torch
from signaturesnet.utilities.io import tensor_to_csv

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utilities.plotting import plot_metric_vs_mutations_classifier_superlow
from utilities.io import csv_to_tensor, read_model

# MIXED DATASET
path = "../../../data/exp_superlow_nummut/classifier/"
num_mut = csv_to_tensor(path + "test_num_mut.csv", device='cpu')
inputs = csv_to_tensor(path + "test_input.csv", device='cpu')
label = csv_to_tensor(path + "test_label.csv", device='cpu').to(torch.int64)

classifier = read_model("../../trained_models/detector")
classifier_guess = classifier(inputs, num_mut)
classification_cutoff = 0.5
classification_results = (classifier_guess >= classification_cutoff).to(torch.int64)
all = torch.cat([num_mut, label, classifier_guess, classification_results], dim=1)
tensor_to_csv(all, '../../../data/exp_superlow_nummut/classifier_results.csv')

plot_metric_vs_mutations_classifier_superlow(classification_results, label, num_mut)


# ROC AUC:
# from sklearn.metrics import roc_curve
# from sklearn.metrics import roc_auc_score

# fpr, tpr, thresh = roc_curve(label.detach().numpy(), classifier_guess.detach().numpy(), pos_label=1)
# print(thresh)
# auc_score = roc_auc_score(label.detach().numpy(), classifier_guess.detach().numpy())
# print(auc_score)

# import matplotlib.pyplot as plt

# # plot roc curves
# plt.plot(fpr, tpr, linestyle='--',color='orange')
# plt.title('ROC curve')
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive rate')
# plt.show()
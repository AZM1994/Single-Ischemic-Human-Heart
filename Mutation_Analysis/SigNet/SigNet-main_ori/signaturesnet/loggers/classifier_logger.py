import os
import sys

import torch
import wandb

from signaturesnet.utilities.metrics import accuracy, false_random, false_realistic


class ClassifierLogger:
    def __init__(self):
        self.threshold = 0.5

    def log(self,
            train_loss,
            train_prediction,
            train_label,
            val_loss,
            val_prediction,
            val_label,
            step):

        wandb.log({"train_loss": train_loss})
        wandb.log({"val_loss": val_loss})

        wandb.log({"train_accuracy": accuracy(train_prediction, train_label).item()})
        wandb.log({"val_accuracy": accuracy(val_prediction, val_label).item()})

        wandb.log({"train_false_realistic": false_realistic(train_prediction, train_label).item()})
        wandb.log({"val_false_realistic": false_realistic(val_prediction, val_label).item()})

        wandb.log({"train_false_random": false_random(train_prediction, train_label).item()})
        wandb.log({"val_false_random": false_random(val_prediction, val_label).item()})

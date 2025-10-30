import torch
import wandb

class NumMutLogger:
    def __init__(self):
        pass

    def log(self,
            train_loss,
            train_classification_metrics,
            val_loss,
            val_classification_metrics,
            step):

        wandb.log({"train_loss": train_loss})
        wandb.log({"val_loss": val_loss})

        for metric in train_classification_metrics:
            wandb.log({"train_" + metric: train_classification_metrics[metric]})

        for metric in val_classification_metrics:
            wandb.log({"val_" + metric: val_classification_metrics[metric]})
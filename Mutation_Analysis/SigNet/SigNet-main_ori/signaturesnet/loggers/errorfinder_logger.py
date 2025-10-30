import torch
import wandb

class ErrorFinderLogger:
    def __init__(self, path, experiment_id):
        pass

    def log(self,
            train_loss,
            pi_metrics_train,
            val_loss,
            pi_metrics_val,
            val_values_lower,
            val_values_upper,
            val_nummut,
            step):

        wandb.log({"train_loss": train_loss})
        wandb.log({"val_loss": val_loss})

        for key in pi_metrics_train.keys():
            wandb.log({"train_%s"%key: pi_metrics_train[key].item()})
        
        for key in pi_metrics_val.keys():
            wandb.log({"val_%s"%key: torch.mean(pi_metrics_val[key]).item()})
import numpy as np
import matplotlib.pyplot as plt
import torch
from torch.utils.tensorboard import SummaryWriter
import wandb

from signaturesnet.utilities.metrics import get_MSE, get_kl_divergence, sets_distances, get_distances_metrics
from signaturesnet.utilities.plotting import get_correlation_matrix


class GeneratorLogger:
    def __init__(self, train_inputs, val_inputs, signatures, device):
        self.counter = 0
        self.plot_freq = 10
        self.metrics = {
            "mse": get_MSE,
        }
        self.train_inputs = train_inputs
        self.presence = torch.mean(train_inputs, dim=0)
        self.val_inputs = val_inputs
        self.signatures = signatures
        self.device = device

    def log(self,
            train_loss,
            train_prediction,
            train_label,
            train_mu,
            train_sigma, 
            val_loss,
            val_prediction,
            val_label,
            val_mu,
            val_sigma,
            model,
            step):

        wandb.log({"train_loss": train_loss})
        wandb.log({"val_loss": val_loss})

        for metric_name in self.metrics.keys():
            metric = self.metrics[metric_name]
            wandb.log({"train_" + metric_name: metric(train_prediction, train_label).item()})
            wandb.log({"val_" + metric_name: metric(val_prediction, val_label).item()})

        # val_mse = get_MSE(val_prediction, val_label).item()

        # Between-examples variance
        train_variance = torch.mean(torch.var(train_prediction, dim=0)).item()
        val_variance = torch.mean(torch.var(val_prediction, dim=0)).item()
        wandb.log({"train_pred_variance": train_variance})
        wandb.log({"val_pred_variance": val_variance})

        # Within-examples variance
        train_variance = torch.mean(torch.var(train_mu, dim=1)).item()
        val_variance = torch.mean(torch.var(val_mu, dim=1)).item()
        wandb.log({"train_mu_variance": train_variance})
        wandb.log({"val_mu_variance": val_variance})
        
        def KL(mu, sigma):
            return (0.5*(sigma + mu**2 - torch.log(sigma) - 1).sum(dim=1)).mean(dim=0)
        
        val_KL = KL(val_mu, val_sigma).item()
        wandb.log({"train_KL": KL(train_mu, train_sigma).item()})
        wandb.log({"val_KL": val_KL})

        wandb.log({"train_mu": torch.mean(train_mu)})
        wandb.log({"val_mu": torch.mean(val_mu)})

        wandb.log({"train_sigma": torch.mean(train_sigma)})
        wandb.log({"val_sigma": torch.mean(val_sigma)})

        train_DQ99R = None
        if self.counter % self.plot_freq == 0:
            n = self.train_inputs.size(0) + self.train_inputs.size(1)
            noise = torch.randn((n, model.latent_dim), device=self.device)
            with torch.no_grad():
                generated = model.decode(noise).detach()
            train_real_dists, train_fake_dists = sets_distances(
                real=self.train_inputs,
                fake=generated[:self.train_inputs.size(0)])
            val_real_dists, val_fake_dists = sets_distances(
                real=self.val_inputs,
                fake=generated[self.train_inputs.size(0):])
            train_real_metrics = get_distances_metrics(train_real_dists)
            train_fake_metrics = get_distances_metrics(train_fake_dists)
            val_real_metrics = get_distances_metrics(val_real_dists)
            val_fake_metrics = get_distances_metrics(val_fake_dists)
            train_DQ99R = train_real_metrics["quantiles"][-1]
            wandb.log({"train_DQ99R" : train_real_metrics["quantiles"][-1],
                       "train_DQ99G" : train_fake_metrics["quantiles"][-1]})
            wandb.log({"val_DQ99R" : val_real_metrics["quantiles"][-1],
                       "val_DQ99G" : val_fake_metrics["quantiles"][-1]})
            plt.switch_backend('agg')
            if self.counter % 50*self.plot_freq == 0:
                fig, _ = get_correlation_matrix(generated, self.signatures.columns[1:])
                wandb.log({"correlation": wandb.Image(fig)})
                plt.cla()
                plt.clf()
                plt.close()

                gen_presence = torch.mean(generated, dim=0)
                presence_errors = (gen_presence - self.presence).detach().cpu().numpy()
                # fig2 = plt.figure()
                plt.bar(np.array(list(range(presence_errors.shape[0]))), presence_errors)
                # plt.show()
                wandb.log({"error_by_signature": plt})
            
        self.counter += 1
        return train_DQ99R
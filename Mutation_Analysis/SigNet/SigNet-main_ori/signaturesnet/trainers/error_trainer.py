
import collections
import copy
import os
import pathlib
import sys

import numpy as np
import pandas as pd

import torch
from torch import linalg
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm
import wandb

from signaturesnet.utilities.data_partitions import DataPartitions
from signaturesnet.utilities.io import save_model, read_model
from signaturesnet.utilities.metrics import distance_to_interval, get_pi_metrics
from signaturesnet.models import ErrorFinder
from signaturesnet.loggers import ErrorFinderLogger
from signaturesnet.modules import CombinedFinetuner
from signaturesnet.utilities.plotting import plot_error_by_sig, plot_width_by_sig

class ErrorTrainer:
    def __init__(self,
                 iterations,
                 train_data,
                 val_data,
                 loging_path="../runs",
                 num_classes=72,
                 model_path=None,
                 device="cuda:0",
                 data_folder="../data"):
        self.iterations = iterations  # Iteration refers to passes through all dataset
        self.num_classes = num_classes
        self.train_dataset = train_data
        self.val_dataset = val_data
        self.device = device
        self.model_path = model_path
        self.logger = ErrorFinderLogger(
            path=None,
            experiment_id=None)
        self.data_folder = data_folder

    def __loss(self,
               label,
               pred_lower,
               pred_upper):
        lagrange_base = float(self.loss_params["lagrange_base"])
        lagrange_high_error_sigs = float(self.loss_params["lagrange_high_error_sigs"])
        lagrange_pnorm = float(self.loss_params["lagrange_pnorm"])
        pnorm_order = int(self.loss_params["pnorm_order"] )
        lagrange_smalltozero = float(self.loss_params["lagrange_smalltozero"])
        pnorm_order = pnorm_order if pnorm_order%2 == 1 else pnorm_order + 1
        
        _EPS = 1e-6
        batch_size = float(pred_lower.shape[0])
        lower = label - pred_lower
        upper = pred_upper - label
        lower = nn.ReLU()(-lower)  # 0 if lower than label
        upper = nn.ReLU()(-upper)  # 0 if higher than label

        interval_length = ((pred_upper - pred_lower)**2)/batch_size
        
        # This is super sketchy
        lagrange_vector = torch.ones(72).to(self.device)*lagrange_base
        lagrange_vector[2] = lagrange_high_error_sigs
        lagrange_vector[4] = lagrange_high_error_sigs
        lagrange_vector[22] = lagrange_high_error_sigs
        lagrange_vector[43] = lagrange_high_error_sigs
        loss_by_mutation_signature =\
            interval_length +\
            lagrange_vector*(lower + upper)

        # p-norm by signature to avoid high errors
        loss_by_mutation = linalg.norm(lagrange_pnorm *\
            loss_by_mutation_signature, ord=pnorm_order, axis=1)
        loss = torch.mean(loss_by_mutation)

        # Send small to 0
        loss += lagrange_smalltozero*\
            torch.mean(torch.abs(pred_upper[label <= _EPS]))
        return loss

    def __meta_loss(self,
                    label,
                    pred_lower,
                    pred_upper):
        _EPS = 1e-6
        batch_size = float(pred_lower.shape[0])
        lower = label - pred_lower
        upper = pred_upper - label
        lower = nn.ReLU()(-lower)  # 0 if lower than label
        upper = nn.ReLU()(-upper)  # 0 if higher than label

        # Manage this case
        # assert(torch.isnan(pred_upper).any() == False)  # int len
        # assert(torch.isnan(pred_lower).any() == False)  # int len

        # Get interval length
        interval_length = torch.mean((pred_upper - pred_lower)**2)

        precision_by_sig = torch.mean((lower + upper < _EPS).to(torch.float32), dim=0)
        high_error_sigs = torch.sum(precision_by_sig<0.95)
        meta_loss = interval_length + high_error_sigs
        return meta_loss

    def objective(self,
                  batch_size,
                  lr,
                  num_neurons_pos,
                  num_neurons_neg,
                  num_hidden_layers_pos,
                  num_hidden_layers_neg,
                  loss_params,
                  plot=False):
        print(batch_size, lr, num_neurons_pos, num_hidden_layers_pos,
              num_neurons_neg, num_hidden_layers_neg, loss_params)
        self.loss_params = loss_params
        dataloader = DataLoader(dataset=self.train_dataset,
                                batch_size=int(batch_size),
                                shuffle=True)
        model = ErrorFinder(num_classes=self.num_classes,
                            num_hidden_layers_pos=int(num_hidden_layers_pos),
                            num_units_pos=int(num_neurons_pos),
                            num_hidden_layers_neg=int(num_hidden_layers_neg),
                            num_units_neg=int(num_neurons_neg))
        model.to(self.device)

        log_freq = 5

        optimizer = optim.Adam(model.parameters(),
                               lr=lr)
        meta_loss_vals = collections.deque(maxlen=20)
        max_found = -np.inf
        step = 0
        for _ in range(self.iterations):
            for _, train_label, train_weight_guess, num_mut, classification in tqdm(dataloader):
                optimizer.zero_grad()
                if "cuda" in str(self.device):
                    train_label = train_label.to(self.device)
                    train_weight_guess = train_weight_guess.to(self.device)
                    num_mut = num_mut.to(self.device)
                    classification = classification.to(self.device)
                train_pred_upper, train_pred_lower = model(weights=train_weight_guess,
                                                           num_mutations=num_mut,
                                                           classification=classification)

                # Compute loss
                train_loss = self.__loss(label=train_label,
                                         pred_lower=train_pred_lower,
                                         pred_upper=train_pred_upper)
                if torch.isnan(train_loss).any():
                    print("NANs appeared while training! Ending training now")
                    return -2.0
                train_loss.backward()
                optimizer.step()

                with torch.no_grad():
                    val_pred_upper, val_pred_lower = model(weights=self.val_dataset.prev_guess,
                                                           num_mutations=self.val_dataset.num_mut,
                                                           classification=self.val_dataset.classification)
                    val_loss = self.__loss(label=self.val_dataset.labels,
                                           pred_lower=val_pred_lower,
                                           pred_upper=val_pred_upper)
                    val_meta_loss = self.__meta_loss(label=self.val_dataset.labels,
                                                     pred_lower=val_pred_lower,
                                                     pred_upper=val_pred_upper)
                    meta_loss_vals.append(val_meta_loss.item())
                    max_found = max(max_found, -np.nanmean(meta_loss_vals))

                if plot and step % log_freq == 0:
                    pi_metrics_train = get_pi_metrics(train_label, train_pred_lower, train_pred_upper)
                    pi_metrics_val = get_pi_metrics(self.val_dataset.labels, val_pred_lower, val_pred_upper, collapse=False, dim=1)
                    self.logger.log(train_loss=train_loss,
                                    pi_metrics_train=pi_metrics_train,
                                    val_loss=val_loss,
                                    pi_metrics_val=pi_metrics_val,
                                    val_values_lower=val_pred_lower,
                                    val_values_upper=val_pred_upper,
                                    val_nummut=self.val_dataset.num_mut,
                                    step=step)
                del train_label, train_weight_guess, num_mut, classification, train_pred_upper, train_pred_lower, val_pred_upper, val_pred_lower
                # if self.model_path is not None and step % 500 == 0:
                #     save_model(model=model, directory=self.model_path)
                step += batch_size
        if self.model_path is not None:
            save_model(model=model, directory=self.model_path)

        with torch.no_grad():
            val_pred_upper, val_pred_lower = model(weights=self.val_dataset.prev_guess,
                                    num_mutations=self.val_dataset.num_mut,
                                    classification=self.val_dataset.classification)
            fig_error = plot_error_by_sig(label=self.val_dataset.labels.cpu(),
                                    pred_upper=val_pred_upper.cpu(),
                                    pred_lower=val_pred_lower.cpu(),
                                    sigs_names=list(pd.read_excel(self.data_folder + "/data.xlsx").columns)[1:])
            wandb.log({"Validation Error by Sig": wandb.Image(fig_error)})

            fig_width = plot_width_by_sig(pred_upper=val_pred_upper.cpu(),
                                        pred_lower=val_pred_lower.cpu(),
                                        sigs_names=list(pd.read_excel(self.data_folder + "/data.xlsx").columns)[1:])
            wandb.log({"Validation Width by Sig": wandb.Image(fig_width)})

        return max_found


def train_errorfinder(config, data_folder="data") -> float:
    from signaturesnet.utilities.io import read_data
    from modules.combined_finetuner import baseline_guess_to_combined_finetuner_guess

    # Select training device
    dev = "cuda" if config["device"] == "cuda" and torch.cuda.is_available() else "cpu"
    device = torch.device(dev)
    print("Using device:", device)

    # Set paths
    errorfinder_path = os.path.join(config["models_dir"], config["model_id"])
    classifier_path = os.path.join(config["models_dir"], config["classifier_id"])
    finetuner_low_path = os.path.join(config["models_dir"], config["finetuner_low_id"])
    finetuner_large_path = os.path.join(config["models_dir"], config["finetuner_large_id"])

    if config["enable_logging"]:
        run = wandb.init(project=config["wandb_project_id"],
                entity='sig-net',
                config=config,
                name=config["model_id"])

    # Load data
    train_real_low, val_real_low = read_data(experiment_id=config["data_id"],
                                            source="low",
                                            data_folder=data_folder,
                                            device=dev)
    train_real_large, val_real_large = read_data(experiment_id=config["data_id"],
                                                source="large",
                                                data_folder=data_folder,
                                                device=dev)
    # Join datasets
    train_data = train_real_low
    train_data.append(train_real_large)
    train_data.perm()

    val_data = val_real_low
    val_data.append(val_real_large)

    model = CombinedFinetuner(low_mum_mut_dir=finetuner_low_path,
                              large_mum_mut_dir=finetuner_large_path,
                              device=dev)
    classifier = read_model(classifier_path, device=dev)

    train_data = baseline_guess_to_combined_finetuner_guess(model=model,
                                                            classifier=classifier,
                                                            data=train_data)
    val_data = baseline_guess_to_combined_finetuner_guess(model=model,
                                                          classifier=classifier,
                                                          data=val_data)

    # train_data.to(device=device)  # We keep the train data in cpu to save memory
    # val_data.to(device=device)  # If still a problem we can make the val_data smaller
    # torch.cuda.empty_cache()

    trainer = ErrorTrainer(iterations=config["iterations"],  # Passes through all dataset
                           train_data=train_data,
                           val_data=val_data,
                           num_classes=config["num_classes"],
                           device=device,
                           model_path=errorfinder_path,
                           data_folder=data_folder)

    min_val = trainer.objective(batch_size=config["batch_size"],
                                lr=config["lr"],
                                num_neurons_pos=config["num_neurons_pos"],
                                num_neurons_neg=config["num_neurons_neg"],
                                num_hidden_layers_pos=config["num_hidden_layers_pos"],
                                num_hidden_layers_neg=config["num_hidden_layers_neg"],
                                loss_params=config["loss_params"],
                                plot=config["enable_logging"])
    if config["enable_logging"]:
        run.finish()
    return min_val
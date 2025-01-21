
import collections
import os
import sys

import numpy as np
import pandas as pd
from pandas.core.algorithms import mode
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from tqdm import tqdm
import wandb

from signaturesnet import DATA
from signaturesnet.utilities.metrics import get_jensen_shannon, get_kl_divergence
from signaturesnet.utilities.io import save_model
from signaturesnet.utilities.generator_data import GeneratorData
from signaturesnet.utilities.oversampler import OverSampler, CancerTypeOverSampler
from signaturesnet.models import Generator
from signaturesnet.loggers.generator_logger import GeneratorLogger

class GeneratorTrainer:
    def __init__(self,
                 iterations,
                 train_data,
                 val_data,
                 signatures,
                 lagrange_param=1.0,
                 loging_path="../runs",
                 num_classes=72,
                 log_freq=100,
                 model_path=None,  # File where to save model learned weights None to not save
                 device=torch.device("cuda:0")):
        self.iterations = iterations  # Now iteration refers to passes through all dataset
        self.num_classes = num_classes
        self.device = device
        self.log_freq = log_freq
        self.lagrange_param = lagrange_param
        self.model_path = model_path
        self.train_dataset = train_data
        self.val_dataset = val_data

        self.logger = GeneratorLogger(
            train_inputs=train_data.inputs,
            val_inputs=val_data.inputs,
            signatures=signatures,
            device=device)

        os = CancerTypeOverSampler(self.train_dataset.inputs, self.train_dataset.cancer_types)
        # os = OverSampler(self.train_dataset.inputs)
        self.train_dataset.inputs = os.get_oversampled_set()

    def __loss(self, input, pred, z_mu, z_std):
        kl_div = (0.5*(z_std.pow(2) + z_mu.pow(2) - 2*torch.log(z_std) - 1).sum(dim=1)).mean(dim=0)
        reconstruction = nn.MSELoss()(input, pred)
        # reconstruction = get_jensen_shannon(predicted_label=pred, true_label=input)
        return reconstruction + self.adapted_lagrange_param*kl_div

    def objective(self,
                  batch_size,
                  lr_encoder,
                  lr_decoder,
                  num_hidden_layers,
                  latent_dim,
                  plot=False):

        print(batch_size, lr_encoder, lr_decoder,
              num_hidden_layers, latent_dim)

        dataloader = DataLoader(dataset=self.train_dataset,
                                batch_size=int(batch_size),
                                shuffle=True)
        model = Generator(input_size=int(self.num_classes),
                          num_hidden_layers=int(num_hidden_layers),
                          latent_dim=int(latent_dim),
                          device=self.device.type)
        model.to(self.device)

        # wandb.watch(model, log_freq=100)

        # optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
        optimizer = optim.Adam([
            {'params': model.encoder_layers.parameters(), 'lr': lr_encoder},
            {'params': model.decoder_layers.parameters(), 'lr': lr_decoder}
        ])

        # l_vals = collections.deque(maxlen=50)
        # max_found = -np.inf
        step = 0
        # total_steps = 1000*len(self.train_dataset)
        total_steps = self.iterations*len(self.train_dataset)
        # self.batch_size_factor = batch_size/len(self.train_dataset)
        self.batch_size_factor = 1.
        train_DQ99R = None
        for iteration in range(self.iterations):
            for train_input in tqdm(dataloader):
                model.train()  # NOTE: Very important! Otherwise we zero the gradient
                optimizer.zero_grad()
                train_pred, train_mean, train_std = model(train_input)
                # self.adapted_lagrange_param = self.lagrange_param
                if step < total_steps*0.8:
                    self.adapted_lagrange_param = self.lagrange_param * \
                        float(total_steps - step)/float(total_steps)
                else:
                    self.adapted_lagrange_param = self.lagrange_param * \
                        float(total_steps - total_steps*0.8)/float(total_steps)
                train_loss = self.__loss(input=train_input,
                                         pred=train_pred,
                                         z_mu=train_mean,
                                         z_std=train_std)

                train_loss.backward()
                optimizer.step()

                model.eval()
                with torch.no_grad():
                    val_pred, val_mean, val_std = model(
                        self.val_dataset.inputs)
                    val_loss = self.__loss(input=self.val_dataset.inputs,
                                           pred=val_pred,
                                           z_mu=val_mean,
                                           z_std=val_std)
                    # l_vals.append(val_loss.item())
                    # max_found = max(max_found, -np.nanmean(l_vals))

                if plot and step % self.log_freq == 0:
                    current_train_DQ99R = self.logger.log(train_loss=train_loss,
                                                        train_prediction=train_pred,
                                                        train_label=train_input,
                                                        val_loss=val_loss,
                                                        val_prediction=val_pred,
                                                        val_label=self.val_dataset.inputs,
                                                        train_mu=train_mean,
                                                        train_sigma=train_std,
                                                        val_mu=val_mean,
                                                        val_sigma=val_std,
                                                        step=step,
                                                        model=model)
                    train_DQ99R = current_train_DQ99R if current_train_DQ99R is not None else train_DQ99R

                if self.model_path is not None and step % 500 == 0:
                    save_model(model=model, directory=self.model_path)
                step += 1
        if self.model_path is not None:
            save_model(model=model, directory=self.model_path)
        
        # Return last mse and KL obtained in validation
        return train_DQ99R

def log_results(config, train_DQ99R, out_csv_path):
    model_results = pd.DataFrame({"batch_size": [config["batch_size"]],
                                  "lr_encoder": [config["lr_encoder"]],
                                  "lr_decoder": [config["lr_decoder"]],
                                  "num_hidden_layers": [config["num_hidden_layers"]],
                                  "latent_dim": [config["latent_dim"]],
                                  "lagrange_param": [config["lagrange_param"]],
                                  "train_DQ99R": [train_DQ99R]})
    model_results.to_csv(out_csv_path,
                         header=False, index=False, mode="a")

def train_generator(config, data_folder=DATA + "/") -> float:
    """Train a classification model and get the validation score

    Args:
        config (dict): Including all the needed args
        to load data, and train the model 
    """
    from signaturesnet.utilities.io import read_data_generator, sort_signatures

    dev = "cuda" if config["device"] == "cuda" and torch.cuda.is_available(
    ) else "cpu"
    print("Using device:", dev)

    if config["enable_logging"]:
        run = wandb.init(project=config["wandb_project_id"],
                   entity='sig-net',
                   config=config,
                   name=config["model_id"])

    train_data, val_data = read_data_generator(device=dev,
                                               data_id=config['data_id'],
                                               cosmic_version=config['cosmic_version'],
                                               data_folder=data_folder,
                                               type=config['type'])

    signatures = sort_signatures(file=data_folder + "data.xlsx",
                                 mutation_type_order=data_folder + "mutation_type_order.xlsx")

    trainer = GeneratorTrainer(iterations=config["iterations"],  # Passes through all dataset
                               train_data=train_data,
                               val_data=val_data,
                               signatures=signatures,
                               lagrange_param=config["lagrange_param"],
                               num_classes=config["num_classes"],
                               device=torch.device(dev),
                               model_path=os.path.join(config["models_dir"], config["model_id"]))

    train_DQ99R = trainer.objective(batch_size=config["batch_size"],
                                    lr_encoder=config["lr_encoder"],
                                    lr_decoder=config["lr_decoder"],
                                    num_hidden_layers=config["num_hidden_layers"],
                                    latent_dim=config["latent_dim"],
                                    plot=config["enable_logging"])

    wandb.log({"train_DQ99R_score": train_DQ99R})

    if config["enable_logging"]:
        run.finish()
    return train_DQ99R

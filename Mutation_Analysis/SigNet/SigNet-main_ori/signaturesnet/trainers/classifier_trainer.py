import collections
import os

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from tqdm import tqdm
import wandb

from signaturesnet.loggers import ClassifierLogger
from signaturesnet.models import Classifier
from signaturesnet.utilities.data_partitions import DataPartitions
from signaturesnet.utilities.io import save_model

class ClassifierTrainer:
    def __init__(self,
                 iterations,
                 train_data,
                 val_data,
                 sigmoid_params = [5000, 2000],
                 loging_path="../runs",
                 num_classes=72,
                 log_freq=100,
                 model_path=None,  # File where to save model learned weights None to not save
                 cutoff=0.5,
                 device=torch.device("cuda:0")):
        self.iterations = iterations  # Now iteration refers to passes through all dataset
        self.num_classes = num_classes
        self.sigmoid_params = sigmoid_params
        self.device = device
        self.log_freq = log_freq
        self.model_path = model_path
        self.train_dataset = train_data
        self.val_dataset = val_data
        self.cutoff = cutoff
        self.logger = ClassifierLogger()

    def __loss(self, prediction, label):
        return nn.BCELoss()(prediction, label)
         
    def objective(self,
                  batch_size,
                  lr,
                  num_hidden_layers,
                  num_neurons,
                  run=None):

        print(batch_size, lr, num_hidden_layers, num_neurons)

        dataloader = DataLoader(dataset=self.train_dataset,
                                batch_size=int(batch_size),
                                shuffle=True)
        model = Classifier(num_hidden_layers=int(num_hidden_layers),
                          num_units=int(num_neurons),
                          sigmoid_params=self.sigmoid_params)
        model.to(self.device)

        optimizer = optim.Adam(model.parameters(), lr=lr)

        l_vals = collections.deque(maxlen=50)
        max_found = -np.inf
        step = 0
        for iteration in range(self.iterations):
            for train_input, train_label, _, num_mut, _ in tqdm(dataloader):
                model.train()  # NOTE: Very important! Otherwise we zero the gradient
                optimizer.zero_grad()                
                train_prediction = model(train_input, num_mut)
                train_loss = self.__loss(prediction=train_prediction,
                                         label=train_label)

                train_loss.backward()
                optimizer.step()

                model.eval()
                with torch.no_grad():
                    val_prediction = model(
                        self.val_dataset.inputs, self.val_dataset.num_mut)
                    val_loss = self.__loss(prediction=val_prediction,
                                        label=self.val_dataset.labels)
                    l_vals.append(val_loss.item())
                    max_found = max(max_found, -np.nanmean(l_vals))

                if run and step % self.log_freq == 0:
                    self.logger.log(train_loss=train_loss,
                                    train_prediction=(train_prediction > self.cutoff).type(torch.int64),
                                    train_label=train_label.type(torch.int64),
                                    val_loss=val_loss,
                                    val_prediction=(val_prediction > self.cutoff).type(torch.int64),
                                    val_label=self.val_dataset.labels.type(torch.int64),
                                    step=step)

                if self.model_path is not None and step % 500 == 0:
                    save_model(model=model, directory=self.model_path)
                step += 1
        if self.model_path is not None:
            save_model(model=model, directory=self.model_path)
        
        if run is not None:
            run.finish()
        return max_found


def train_classifier(config) -> float:
    """Train a classification model and get the validation score

    Args:
        config (dict): Including all the needed args
        to load data, and train the model 
    """
    from signaturesnet.utilities.io import read_data_classifier

    dev = "cuda" if config["device"] == "cuda" and torch.cuda.is_available() else "cpu"
    print("Using device:", dev)

    run = None
    if config["enable_logging"]:
        run = wandb.init(project=config["wandb_project_id"],
                        entity='sig-net',
                        config=config,
                        name=config["model_id"])

    train_data, val_data = read_data_classifier(device=dev,
                                                experiment_id=config["data_id"])

    trainer = ClassifierTrainer(iterations=config["iterations"],  # Passes through all dataset
                               train_data=train_data,
                               val_data=val_data,
                               sigmoid_params=config["sigmoid_params"],
                               device=torch.device(dev),
                               model_path=os.path.join(config["models_dir"], config["model_id"]))

    min_val = trainer.objective(batch_size=config["batch_size"],
                                lr=config["lr"],
                                num_hidden_layers=config["num_hidden_layers"],
                                num_neurons=config["num_neurons"],
                                run=run)

    return min_val
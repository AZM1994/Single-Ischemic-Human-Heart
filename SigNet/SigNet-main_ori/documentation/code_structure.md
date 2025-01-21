# Source code structure

- `signaturesnet/`: Contains all src scripts needed to run and train the SigNet models.
  - `configs/`: Contains the model architecture and training parameters of each neural network used in SigNet.
  - `data/`: Refer [here](data_folder.md).
  - `generate_data/`: Collection of scripts used to create the training sets for each module.
  - `hyperparam_optimizers/`: Collection of scripts to find the best network training hyperparameters for each module.
  - `loggers/`: Loggers used to store into [wandb](https://wandb.ai/site) training for each module.
  - `models/`: Holds the neural network architecture definition for each module.
  - `modules/`: Collection of common ANN building blocks.
  - `trained_models/`: Weights of the pre-trained models.
  - `trainers/`: Main training loop of each of SigNet's modules.
  - `train_*.py`: Code to train each particular module.

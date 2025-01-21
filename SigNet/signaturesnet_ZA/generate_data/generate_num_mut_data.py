
import numpy as np
import pandas as pd
import pathlib

from signaturesnet import DATA

def shuffle(inputs, labels):
    indexes = np.random.permutation(inputs.shape[0])
    return inputs.iloc[indexes, :], labels[indexes]

def encoding(muts):
    muts_new = muts.copy()
    muts_new[(muts>=1) & (muts<1.5)] = 0
    muts_new[(muts>=1.5) & (muts<2)] = 1
    muts_new[(muts>=2) & (muts<2.5)] = 2
    muts_new[(muts>=2.5) & (muts<3)] = 3
    muts_new[(muts>=3) & (muts<3.5)] = 4
    muts_new[(muts>=3.5) & (muts<4)] = 5
    muts_new[(muts>=4) & (muts<4.5)] = 6
    muts_new[(muts>=4.5) & (muts<5)] = 7
    muts_new[muts>5] = 8
    return muts_new

if __name__ == "__main__":

    # Read
    real_data_weights = pd.read_csv(DATA + "/real_data/sigprofiler_not_norm_PCAWG.csv", header=0, index_col=0).reset_index(drop=True)
    inputs = pd.concat([real_data_weights, pd.DataFrame(np.zeros((real_data_weights.shape[0], 7)))], axis=1, ignore_index=True)

    real_data = pd.read_csv(DATA + '/real_data/PCAWG_sigProfiler_SBS_signatures_in_samples_v3.csv', sep=',')
    real_data = real_data.iloc[:,3:]
    total_muts = np.log10(real_data.sum(axis=1))

    labels = encoding(total_muts)
    
    import matplotlib.pyplot as plt
    plt.hist(total_muts, bins=10)
    plt.show()
    plt.hist(labels, bins=20)
    plt.show()

    inputs, labels = shuffle(inputs, labels)
    
    # Partition datasets
    n = inputs.shape[0]
    train_inputs = inputs.iloc[:int(n*0.7), :]
    train_labels = labels.iloc[:int(n*0.7)]

    val_inputs = inputs.iloc[int(n*0.7):int(n*0.90), :]
    val_labels = labels.iloc[int(n*0.7):int(n*0.90)]

    test_inputs = inputs.iloc[int(n*0.90):, :]
    test_labels = labels.iloc[int(n*0.90):,]

    # store everything
    pathlib.Path(DATA + '/datasets/num_muts/').mkdir(parents=True, exist_ok=True)
    train_inputs.to_csv(DATA + '/datasets/num_muts/train_input.csv', header=False, index=False)
    train_labels.to_csv(DATA + '/datasets/num_muts/train_label.csv', header=False, index=False)

    val_inputs.to_csv(DATA + '/datasets/num_muts/val_input.csv', header=False, index=False)
    val_labels.to_csv(DATA + '/datasets/num_muts/val_label.csv', header=False, index=False)

    test_inputs.to_csv(DATA + '/datasets/num_muts/test_input.csv', header=False, index=False)
    test_labels.to_csv(DATA + '/datasets/num_muts/test_label.csv', header=False, index=False)

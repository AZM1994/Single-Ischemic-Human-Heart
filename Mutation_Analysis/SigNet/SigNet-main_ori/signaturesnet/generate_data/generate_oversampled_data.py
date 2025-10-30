from signaturesnet.utilities.io import read_data_generator, read_signatures, tensor_to_csv
from signaturesnet.utilities.oversampler import CancerTypeOverSampler
from signaturesnet.utilities.data_generator import DataGenerator
from signaturesnet.utilities.data_partitions import DataPartitions
from signaturesnet.models import Baseline
from signaturesnet import DATA

if __name__ == "__main__":
    train_data, val_data = read_data_generator(device="cpu",
                                               data_id="/real_data",
                                               data_folder=DATA,
                                               cosmic_version='v3',
                                               type='real',
                                               prop_train = 0.8)
    oversampler = CancerTypeOverSampler(train_data.inputs, train_data.cancer_types)
    train_label = oversampler.get_even_set()         # Oversample to create set with same number of samples per cancer type
    val_label = val_data.inputs

    # Create inputs associated to the labels
    signatures = read_signatures(DATA + "/data.xlsx", mutation_type_order=DATA + "/mutation_type_order.xlsx")
    data_generator = DataGenerator(signatures=signatures,
                                   seed=0,
                                   shuffle=True)
    train_input, train_label = data_generator.make_input(labels=train_label,
                                                         split="train",
                                                         large_low="superlow")
    val_input, val_label = data_generator.make_input(labels=val_label,
                                                     split="valid",
                                                     large_low="superlow")
    
    # Run Baseline
    sf = Baseline(signatures)
    train_baseline = sf.get_weights_batch(train_input, n_workers=2)
    val_baseline = sf.get_weights_batch(val_input, n_workers=2)
    
    # Create DataPartitions
    train_data = DataPartitions(inputs=train_input,
                                prev_guess=train_baseline,
                                labels=train_label)
    val_data = DataPartitions(inputs=val_input,
                                prev_guess=val_baseline,
                                labels=val_label)

    tensor_to_csv(train_data.inputs, "../data/exp_superlow_nummut/refitter/train_%s_input.csv"%"superlow")
    tensor_to_csv(train_label, "../data/exp_superlow_nummut/refitter/train_%s_label.csv"%"superlow")
    tensor_to_csv(train_data.prev_guess, "../data/exp_superlow_nummut/refitter/train_%s_baseline.csv"%"superlow")

    tensor_to_csv(val_data.inputs, "../data/exp_superlow_nummut/refitter/val_%s_input.csv"%"superlow")
    tensor_to_csv(val_label, "../data/exp_superlow_nummut/refitter/val_%s_label.csv"%"superlow")
    tensor_to_csv(val_data.prev_guess, "../data/exp_superlow_nummut/refitter/val_%s_baseline.csv"%"superlow")
        
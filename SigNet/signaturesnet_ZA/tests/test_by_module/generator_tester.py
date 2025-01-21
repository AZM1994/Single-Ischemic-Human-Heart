import json
import os
import sys

import pandas as pd
import torch
from pprint import pprint

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utilities.plotting import plot_correlation_matrix, plot_histograms
from utilities.metrics import prop_distances, sets_distances, get_distances_metrics
from utilities.io import csv_to_tensor, read_data_generator, read_model, read_signatures, sort_signatures
from utilities.oversampler import CancerTypeOverSampler

def read_real_data():
    inputs = data_folder + "real_data/PCAWG_data.csv"
    labels = data_folder + "real_data/sigprofiler_not_norm_PCAWG.csv"
    inputs = csv_to_tensor(file=inputs, header=0, index_col=0)
    labels = csv_to_tensor(labels, header=0, index_col=0)
    labels = labels/torch.sum(labels, axis=1).reshape(-1, 1)
    labels = torch.cat([labels, torch.zeros(labels.size(0), 7).to(labels)], dim=1)

    nummut = torch.sum(inputs, dim=1)
    # inputs = normalize_data(inputs,
                            # opportunities_name_or_path="../../data/real_data/3mer_WG_hg37.txt")
                            # opportunities_name_or_path="../../data/real_data/abundances_trinucleotides.txt")
                            # opportunities_name_or_path="../../data/real_data/norm_38.txt")
                            # opportunities_name_or_path="../../data/real_data/new_norm.txt")
    inputs = inputs/torch.sum(inputs, axis=1).view(-1, 1)
    return inputs, labels, nummut

if __name__=="__main__":
    data_folder = "../../../data/"

    real_inputs, real_labels, real_nummut = read_real_data()

    results_dict = {'model_id': [],
                    'DQ99G': [],
                    'DMeanG': [],
                    'DQ99R': [],
                    'DMeanR': [],
                    'prop_distance': [],
                    'mean_prop_distance': []}

    # list_of_dirs = ["../../../trained_models/fix_generator_oversample/cancer_type_oversampler_" + str(i+1) for i in range(20)]
    list_of_dirs = ["../../../trained_models/fix_generator_oversample/fix_generator_light_oversample_bayesian"]
    for directory in list_of_dirs:
        generator_vae = read_model(directory)

        init_args_file = os.path.join(directory, 'init_args.json')
        with open(init_args_file, 'r') as fp:
            init_args = json.load(fp)
        latent_dim = init_args["latent_dim"]

        noise = torch.randn((2000, latent_dim))
        synt_labels_vae = generator_vae.decode(noise)

        real_dists, fake_dists = sets_distances(real_labels, synt_labels_vae)

        real_dists_metrics = get_distances_metrics(real_dists)
        pprint(real_dists_metrics)
        
        fake_dists_metrics = get_distances_metrics(fake_dists)
        pprint(fake_dists_metrics)

        results_dict['model_id'].append(directory.split('/')[-1])
        results_dict['DQ99G'].append(fake_dists_metrics['quantiles'][0][0])
        results_dict['DMeanG'].append(fake_dists_metrics['mean'])
        results_dict['DQ99R'].append(real_dists_metrics['quantiles'][0][0])
        results_dict['DMeanR'].append(real_dists_metrics['mean'])
        results_dict['prop_distance'].append(prop_distances(real_labels, synt_labels_vae)[0])
        results_dict['mean_prop_distance'].append(prop_distances(real_labels, synt_labels_vae)[1])

    results_df = pd.DataFrame(results_dict)
    print(results_df)
    # results_df.to_csv("generator_performance.csv")
    data_dict = {
        "fake": fake_dists,
        "real": real_dists
    }
    plot_histograms(data_dict)

    print(prop_distances(real_labels, synt_labels_vae))

    # Correlation matrices
    signatures = sort_signatures(file=data_folder + "data.xlsx",
                                 mutation_type_order=data_folder + "mutation_type_order.xlsx")
    plot_correlation_matrix(data=real_labels, signatures=signatures)
    plot_correlation_matrix(data=synt_labels_vae, signatures=signatures)

    train_data, val_data = read_data_generator(device='cpu',
                                               data_id="real_data",
                                               cosmic_version='v3',
                                               data_folder=data_folder,
                                               type="real")

    os = CancerTypeOverSampler(train_data.inputs, train_data.cancer_types)
    real_labels_light_augmented = os.get_oversampled_set()
    real_labels_augmented = os.get_even_set()

    plot_correlation_matrix(data=real_labels_light_augmented, signatures=signatures)
    plot_correlation_matrix(data=real_labels_augmented, signatures=signatures)
    # plot_correlation_matrix(data=synt_labels_gan, signatures=signatures)


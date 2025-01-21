import pandas as pd
from signaturesnet import DATA
from signaturesnet.modules.signet_module import SigNet

## for mutation signature analysis
# PTA_Cases
mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/Mutation_Analysis/Hypoxia_PTA_Cases/sSNV_figures/mut_mat_estimiated.csv", header=0, index_col=0)
signet_mut = SigNet(opportunities_name_or_path = "genome")
# results_mut = signet_mut(mutation_dataset = mutations_mut, cutoff = 0.01, only_NNLS = False)
results_mut = signet_mut(mutation_dataset = mutations_mut)
w_mut, u_mut, l_mut, c_mut, _ = results_mut.get_output()
# results_mut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/PTA_Cases_0812')
results_mut.save(path='/Users/zhemingan/Documents/BCH_research/Mutation_Analysis/Hypoxia_PTA_Cases/sSNV_figures/SigNet/signature_weights')

# PTA_IPSC
mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/Mutation_Analysis/Hypoxia_PTA_IPSC/sSNV_figures/mut_mat_estimiated.csv", header=0, index_col=0)
signet_mut = SigNet(opportunities_name_or_path = "genome")
# results_mut = signet_mut(mutation_dataset = mutations_mut, cutoff = 0.01, only_NNLS = False)
results_mut = signet_mut(mutation_dataset = mutations_mut)
w_mut, u_mut, l_mut, c_mut, _ = results_mut.get_output()
# results_mut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/PTA_Cases_0812')
results_mut.save(path='/Users/zhemingan/Documents/BCH_research/Mutation_Analysis/Hypoxia_PTA_IPSC/sSNV_figures/SigNet/signature_weights')



## for enrichment analysis

mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_modified_est.csv", header=0, index_col=0)
# mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_modified.csv", header=0, index_col=0)
signet_mut = SigNet(opportunities_name_or_path = "genome")
# results_mut = signet_mut(mutation_dataset = mutations_mut, cutoff = 0.01, only_NNLS = False)
results_mut = signet_mut(mutation_dataset = mutations_mut)
w_mut, u_mut, l_mut, c_mut, _ = results_mut.get_output()
results_mut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/mutation_mut_mat_all_by_case_decile_est')

mutations_permut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_modified_est.csv", header=0, index_col=0)
# mutations_permut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_modified.csv", header=0, index_col=0)
signet_permut = SigNet(opportunities_name_or_path = "genome")
results_permut = signet_permut(mutation_dataset = mutations_permut)
w_permut, u_permut, l_permut, c_permut, _ = results_permut.get_output()
results_permut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/permutation_mut_mat_all_by_case_decile_est')

## decomposite denovo signatures to COSMIC
mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/SigNet/denovo_sigs.csv", header=0, index_col=0)
signet_mut = SigNet(opportunities_name_or_path = "genome")
results_mut = signet_mut(mutation_dataset = mutations_mut)
w_mut, u_mut, l_mut, c_mut, _ = results_mut.get_output()
results_mut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/denovo_sigs')

## decomposite mut_mat to denovo signatures
mutations_mut = pd.read_csv("/Users/zhemingan/Documents/BCH_research/Hypoxia_Project_Integration/Mutation_Analysis/Hypoxia_PTA_Cases/sSNV_figures/mut_mat_estimiated_PTA_Cases.csv", header=0, index_col=0)
signet_mut = SigNet(opportunities_name_or_path = "genome")
results_mut = signet_mut(mutation_dataset = mutations_mut)
w_mut, u_mut, l_mut, c_mut, _ = results_mut.get_output()
results_mut.save(path='/Users/zhemingan/Documents/BCH_research/SigNet/denovo_sigs_mut_mat')

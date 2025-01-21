import pandas as pd
from signaturesnet import DATA
from signaturesnet.modules.signet_module import SigNet

# Read your mutational data
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/PTA_cases/mut_mat_raw.csv", header=0, index_col=0)
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/PTA_cases/mut_mat_raw_conditional.csv", header=0, index_col=0)
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/PTA_cases/mut_mat_SigNet.csv", header=0, index_col=0)
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/PTA_cases/mut_mat_SigNet_conditional.csv", header=0, index_col=0)

mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/mutation_mut_mat_all_by_case_decile.csv", header=0, index_col=0)
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/permutation_mut_mat_all_by_case_decile.csv", header=0, index_col=0)
# mutations = pd.read_csv("/Users/zhemingan/Documents/SigNet/permutation_mut_mat_all_modified.csv", header=0, index_col=0)
# mutations = mutations.iloc[49:50]

# Load & Run signet
signet = SigNet(opportunities_name_or_path = "genome")
results = signet(mutation_dataset = mutations, cutoff = 0.01, only_NNLS = False)

# Extract results
w, u, l, c, _ = results.get_output()

# Store results
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/retrained/raw')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/retrained/raw_conditional')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/retrained/SigNet')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/retrained/SigNet_conditional')

# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/original/raw')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/original/raw_conditional')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/original/SigNet')
# results.save(path='/Users/zhemingan/Documents/SigNet/PTA_cases/original/SigNet_conditional')

results.save(path='/Users/zhemingan/Documents/SigNet/mutation_mut_mat_all_by_case_decile')
# results.save(path='/Users/zhemingan/Documents/SigNet/permutation_mut_mat_all_by_case_decile')

# Plot figures
# results.plot_results(compute=True, save=True, path='/Users/zhemingan/Documents/SigNet/test')

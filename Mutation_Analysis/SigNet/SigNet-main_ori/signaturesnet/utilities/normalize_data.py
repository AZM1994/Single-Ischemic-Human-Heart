import numpy as np
import pandas as pd
import torch

from signaturesnet import DATA

def normalize_data(input_file, opportunities_name_or_path):
    if opportunities_name_or_path == 'genome':
        return input_file
    else:
        genome_abundances = create_opportunities(DATA + '/real_data/3mer_WG_hg37.txt')
        genome_abundances_tensor = torch.Tensor(genome_abundances)

        if opportunities_name_or_path == 'exome':
            abundances = create_opportunities(DATA + '/real_data/abundances_trinucleotides.txt')
        elif opportunities_name_or_path == 'genome_hg38':
            abundances = create_opportunities(DATA + '/real_data/3mer_WG_hg38.txt')
        else:
            abundances = create_opportunities(opportunities_name_or_path)

        abundances_tensor = torch.Tensor(abundances)
        abundances_tensor = (abundances_tensor/genome_abundances_tensor)/torch.sum(abundances_tensor/genome_abundances_tensor)
        return torch.div(input_file, abundances_tensor)



def create_opportunities(abundances_path):
    # Create opportunities
    order = ['ACA','ACC','ACG','ACT','CCA','CCC','CCG','CCT','GCA','GCC','GCG','GCT','TCA','TCC','TCG','TCT',
             'ACA','ACC','ACG','ACT','CCA','CCC','CCG','CCT','GCA','GCC','GCG','GCT','TCA','TCC','TCG','TCT',
             'ACA','ACC','ACG','ACT','CCA','CCC','CCG','CCT','GCA','GCC','GCG','GCT','TCA','TCC','TCG','TCT',
             'ATA','ATC','ATG','ATT','CTA','CTC','CTG','CTT','GTA','GTC','GTG','GTT','TTA','TTC','TTG','TTT',
             'ATA','ATC','ATG','ATT','CTA','CTC','CTG','CTT','GTA','GTC','GTG','GTT','TTA','TTC','TTG','TTT',
             'ATA','ATC','ATG','ATT','CTA','CTC','CTG','CTT','GTA','GTC','GTG','GTT','TTA','TTC','TTG','TTT']
    
    opp_file = open(abundances_path, 'r')
    opp_file_lines = opp_file.readlines()
    opp_dic = {}
    for line in opp_file_lines:
        try:
            opp_dic[line.split("\t")[0]] = float(line.strip("\n").split("\t")[1])
        except:
            opp_dic[line.split(" ")[0]] = float(line.strip("\n").split(" ")[1])
    opp = np.zeros(96)
    for i in range(len(order)):
        tri = order[i]
        opp[i] = opp_dic[tri]
        b1 = tri[0]
        b2 = tri[1]
        b3 = tri[2]
        complement_base = complement(b3) + complement(b2) + complement(b1)
        opp[i] = opp[i] + opp_dic[complement_base]
    return opp



def complement(base):
    if base == "A":
        return "T"
    if base == "C":
        return "G"
    if base == "G":
        return "C"
    if base == "T":
        return "A"

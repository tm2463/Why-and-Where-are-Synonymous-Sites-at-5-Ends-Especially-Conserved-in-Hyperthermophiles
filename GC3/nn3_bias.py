import numpy as np
import pandas as pd
from collections import defaultdict

start = input('enter concatenated start.csv (codon usage): ').strip()
core = input('enter concatenated core.csv (codon usage): ').strip()

codon_dict = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def make_dicts(region):
    region = pd.read_csv(region)
    meso_dict = dict(zip(region['codon'], region['mesophiles']))
    thermo_dict = dict(zip(region['codon'], region['thermophiles']))
    hyper_dict = dict(zip(region['codon'], region['hyperthermophiles']))
    return meso_dict, thermo_dict, hyper_dict


meso_start_dict, thermo_start_dict, hyper_start_dict = make_dicts(start)
meso_core_dict, thermo_core_dict, hyper_core_dict = make_dicts(core)


def nn3_bias(codon_dict):
    nt_totals = defaultdict(list)
    for codon, value in codon_dict.items():
        third_base = codon[2]
        nt_totals[third_base].append(value)

    nt_bias = {}
    for nt, values in nt_totals.items():
        nt_bias[nt] = np.mean(nt_totals[nt])

    return nt_bias


#create dicts containing nn3 bias
meso_start = nn3_bias(meso_start_dict)
thermo_start = nn3_bias(thermo_start_dict)
hyper_start = nn3_bias(hyper_start_dict)

meso_core = nn3_bias(meso_core_dict)
thermo_core = nn3_bias(thermo_core_dict)
hyper_core = nn3_bias(hyper_core_dict)

plot_data = {"Mesophile 5'": meso_start, "Thermophile 5'": thermo_start, "Hyperthermophile 5'": hyper_start,
             'Mesophile Core': meso_core, 'Thermophile Core': thermo_core, 'Hyperthermophile Core': hyper_core}

nn3_df = pd.DataFrame(plot_data)
nn3_df = pd.DataFrame(nn3_df).T[['A', 'T', 'C', 'G']] #transpose for graph

nn3_df.to_csv('nn3.csv')

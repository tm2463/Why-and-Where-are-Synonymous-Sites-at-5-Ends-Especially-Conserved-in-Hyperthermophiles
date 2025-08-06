#manually combine output csv's in excel
import os
from collections import defaultdict
import pandas as pd

growth_file = input('enter growth file path: ').strip()
cds_dir = input('enter cds directory path: ').strip()
pwd = os.getcwd()
output_dir = os.path.join(pwd, 'codon_usage_results')
os.makedirs('codon_usage_results', exist_ok=True)

#load growth temps file
with open(growth_file) as f:
    growth_temps = f.read().split('\n')
    growth_dict = {}

    for x in growth_temps:
        if x[:7] == 'timeout':
            continue

        try:
            y = x.split('\t')
            genus = y[0]
            temp = y[1]
            growth_dict[genus] = float(temp)
        except:
            print(f'index error at {x}')

cds_list = os.listdir(cds_dir)
try:
    cds_list.remove('.DS_Store')
except:
    print('no ".DS_Store"')

mesophiles, thermophiles, hyperthermophiles = [], [], []
for x in cds_list:
    try:
        if growth_dict[x[:-4]] < 50:
            mesophiles.append(x)
        elif growth_dict[x[:-4]] > 50 and growth_dict[x[:-4]] < 80:
            thermophiles.append(x)
        elif growth_dict[x[:-4]] > 80:
            hyperthermophiles.append(x)
    except:
        continue

cds_list = []
cds_list.append(mesophiles)
cds_list.append(thermophiles)
cds_list.append(hyperthermophiles)

mesophile_cds, thermophile_cds, hyperthermophile_cds = [], [], []
for group, new_list in zip(cds_list, [mesophile_cds, thermophile_cds, hyperthermophile_cds]):
    for genus in group:
        with open(os.path.join(cds_dir, f'{genus}')) as f:
            data = f.read().split('>')
            cds = []

            for x in data:
                seq = ''.join(x.split('\n')[1:])
                cds.append(seq)

        new_list.append(cds)

cds = []
cds.append(mesophile_cds)
cds.append(thermophile_cds)
cds.append(hyperthermophile_cds)

#collect starting 10 codons and middle 10 codons for each gene in each temp group
starts, cores = [], []
for group in cds:
    group_start = []
    group_core = []

    for genera in group:
        start, core = [], []

        for gene in genera:
            try:
                codon_length = len(gene)/3

                if gene[3:33]:
                    start.append(gene[3:33])

                mid = codon_length//2
                core_start = int((mid-5)*3)
                core_end = int((mid+5)*3)
                if gene[core_start:core_end]:
                    core.append(gene[core_start:core_end])

            except:
                continue

        group_start.append(start)
        group_core.append(core)

    starts.append(group_start)
    cores.append(group_core)

codons = {
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

synonymous_counts = {
    'F': 2, 'L': 6, 'S': 6, 'Y': 2,
    'C': 2, 'W': 1, 'P': 4, 'H': 2,
    'Q': 2, 'R': 6, 'I': 3, 'M': 1,
    'T': 4, 'N': 2, 'K': 2, 'V': 4,
    'A': 4, 'D': 2, 'E': 2, 'G': 4
}


#write csv files for counting amino acids and codon usage for 5' ends and gene cores
def make_csv(region, filename):
    for group, names in zip(region, ['meso', 'thermo', 'hyper']):
        codon_count = defaultdict(int)
        aa_count = defaultdict(int)

        for genera in group:
            for gene in genera:
                for x in range(0, len(gene), 3):
                    codon_count[gene[x:x+3]] += 1
                    aa_count[codons[gene[x:x+3]]] += 1

        df = pd.DataFrame(codon_count.items())
        df = df.sort_values(by=0)

        codon = [codons[aa] for aa in df[0]]

        df[2] = codon

        #ignore stop, methionine and tryptophan codons in codon analysis
        df = df.iloc[:, [2, 0, 1]]
        df = df[df[2] != '*']
        df = df[df[2] != 'M']
        df = df[df[2] != 'W']

        df = df.sort_values(by=2) #order amino acids alphabetically for easy 'eyeball' comparisons

        #ignore stop, methionine and tryptophan codons in amino acid analysis
        aa_df = pd.DataFrame(aa_count.items())
        aa_df = aa_df[aa_df[0] != '*']
        aa_df = aa_df[aa_df[0] != 'M']
        aa_df = aa_df[aa_df[0] != 'W']

        #write csvs
        df.to_csv(os.path.join(output_dir, f'{names}_{filename}.csv'), index=False)
        aa_df.to_csv(os.path.join(output_dir, f'{names}_{filename}_count.csv'), index=False)

make_csv(starts, 'start')
make_csv(cores, 'core')

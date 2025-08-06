import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from tqdm import tqdm

permutations = 1000 #edit value to change no. permutations

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

aa_dict = defaultdict(list)
for codon, aa in codon_dict.items():
    aa_dict[aa].append(codon)


#count intergenic nucleotide frequencies
def count_nucleotides(seq):
    a_count = seq.count('A')
    t_count = seq.count('T')
    c_count = seq.count('C')
    g_count = seq.count('G')

    total = a_count + t_count + c_count + g_count

    a_freq = a_count / total
    t_freq = t_count / total
    c_freq = c_count / total
    g_freq = g_count / total
    return (a_freq, t_freq, c_freq, g_freq)


#create utr permutations
def utr_permutations(freqs, num=int(permutations)):
    nucleotides = np.array(['A', 'T', 'C', 'G'])
    probs = np.array(freqs)
    #create array and fill with randomly nucleotides weighted by intergenic
    utr = np.random.choice(nucleotides, size=(num, 30), p=probs)
    output = []
    for row in utr:
        seqs = ''.join(row)
        output.append(seqs)
    return output


#create 5' synonymous mutations
def synonymous_permutation(seq, freqs, num=int(permutations)):
    out = np.empty((num, 20), dtype='U3') #3 character strings for codons

    probs = {}
    probs['A'] = freqs[0]
    probs['T'] = freqs[1]
    probs['C'] = freqs[2]
    probs['G'] = freqs[3]

    for x in range(20):
        codon = seq[x*3:(x+1)*3]
        aa = codon_dict.get(codon) #call alternative codons
        synonymous = aa_dict[aa]

        if len(synonymous) == 1: #proceed to next codon in case of M or W
            out[:, x] = synonymous[0]
            continue

        #bias synonymous codon choice based on intergenic nt freqs
        weights = []
        for codon in synonymous:
            weighting = probs[codon[0]] * probs[codon[1]] * probs[codon[2]]
            weights.append(weighting)

        weights = np.array(weights)
        total_weight = weights.sum()
        weights = weights/total_weight #normalise

        out[:, x] = np.random.choice(synonymous, size=num, p=weights) #store weighted permutations

    #format permutations for later use
    output = []
    for row in out:
        seqs = ''.join(row)
        output.append(seqs)
    return output


def process_gene(n, gene, freqs, permutations, gen_perm_dir):
    try:
        lines = gene.split('\n')
        seq = lines[1][30:] #split 5' UTR and 5' CDS
    except:
        return

    utr_perms = utr_permutations(freqs, permutations)
    cds_perms = synonymous_permutation(seq, freqs, permutations)

    with open(os.path.join(gen_perm_dir, f'{n}.txt'), 'w') as out_file:
        for n, (utr, cds) in enumerate(zip(utr_perms, cds_perms)):
            string = utr + cds
            out_file.write(f'>{n}\n{string}\n')


def main():
    pwd = input('enter elongated sequences path: ').strip()
    inter = input('enter intergenome path: ').strip()
    cwd = os.getcwd()
    os.makedirs(os.path.join(cwd, 'shuffle_output'), exist_ok=True)
    out_dir = os.path.join(cwd, 'shuffle_output')

    inter_dict = defaultdict(list)
    inter_list = os.listdir(inter)
    try:
        inter_list.remove('.DS_Store')
    except:
        pass

    for x in inter_list:
        with open(os.path.join(inter, x)) as f:
            file = f.read().strip()
            inter_dict[x[:-4]] += count_nucleotides(file)

    genome_list = os.listdir(f'{pwd}')

    for genome in tqdm(genome_list, desc='genera completed'):
        gen_perm_dir = os.path.join(out_dir, f'{genome[:-4]}')
        os.makedirs(gen_perm_dir, exist_ok=True)

        freqs = inter_dict[genome[:-4]]

        with open(os.path.join(pwd, genome)) as f:
            genes = f.read().split('>')[1:]

        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            for n, gene in enumerate(genes):
                executor.submit(process_gene, n, gene, freqs, permutations, gen_perm_dir)
                

if __name__ == '__main__':
    main()
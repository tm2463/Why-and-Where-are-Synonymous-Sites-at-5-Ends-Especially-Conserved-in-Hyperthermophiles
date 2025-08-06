import os
import pandas as pd
from scipy.stats import ttest_ind

input_dir = input('enter rnaseq master directory: ').strip()
growth_file = input('enter growth file path: ').strip()
pwd = os.getcwd()

input_list = os.listdir(input_dir)
if '.DS_Store' in input_list:
    input_list.remove('.DS_Store')

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


def get_gc(seq):
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


with (open(os.path.join(pwd, 'output.csv'), 'w') as output):
    output.write('genus,ogt,sig,p\n')

    for x in input_list:
        genera = x.split('_')[0]
        print(f'{x} -> OGT: {growth_dict[genera]}°C')

        with open(os.path.join(input_dir, x, 'salmon_output', 'quant.sf')) as f:
            quant_df = pd.read_csv(f, sep='\t')
            tpm_dict = dict(zip(quant_df['Name'], quant_df['TPM']))

        gc_dict = {}
        reference_list = os.listdir(os.path.join(input_dir, x))
        reference = ''.join([y for y in reference_list if y.endswith('.fna')])

        with open(os.path.join(input_dir, x, reference)) as r:
            ref_file = r.read().strip().split('>')

        for gene in ref_file:
            data = gene.split('\n')

            try:
                header = data[0].split(' ')[0].strip()
                cds = data[1][2:32:3]
                gc_dict[header] = get_gc(cds)
            except:
                continue

        sorted_genes = sorted(tpm_dict.items(), key=lambda item: item[1])
        cutoff = int(len(sorted_genes) * 0.05)
        top_tpm = dict(sorted_genes[-cutoff:])
        bottom_tpm = dict(sorted_genes[:cutoff])

        top_gc = []
        for top in top_tpm:
            if top in gc_dict:
                top_gc.append(gc_dict[top])
            else:
                print(f'missing {top}')

        bottom_gc = []
        for bottom in bottom_tpm:
            if bottom in gc_dict:
                bottom_gc.append(gc_dict[bottom])
            else:
                print(f'missing {bottom}')

        t, p = ttest_ind(top_gc, bottom_gc)

        if p < 0.05:
            sig = 1
        else:
            sig = 0

        print(f"t = {t}, p = {p}")
        output.write(f'{x},{growth_dict[genera]},{sig},{p}\n')

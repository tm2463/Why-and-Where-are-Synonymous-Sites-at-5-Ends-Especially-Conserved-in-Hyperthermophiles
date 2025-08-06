import os
import numpy as np
from Bio import SeqIO


def get_gc_content(seq):
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def string_parse(string):
    s = string.strip()
    s = s.strip('[]')

    string = []
    for x in s.split(','):
        try:
            string.append(float(x.strip()))
        except:
            continue

    if len(string) == 5:
        return string
    else:
        return


directory = os.getcwd()
growth_file = input('enter predicted OGT file path: ').strip()

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

genome_dir = input('enter genome directory path: ').strip()
genome_list = os.listdir(genome_dir)

if '.DS_Store' in genome_list:
    genome_list.remove('.DS_Store')

gc_contents = {}
for x, genome in enumerate(genome_list):
        try:
            with open(os.path.join(genome_dir, genome)) as f:
                genome_seq = []

                for record in SeqIO.parse(f, 'fasta'):
                    genome_seq.append(str(record.seq))

                try:
                    if growth_dict[genome[:-4]]:
                        genome_seq = ''.join(genome_seq)
                        gc_content = get_gc_content(genome_seq)
                        gc_contents[genome[:-4]] = gc_content

                except:
                    continue

        except:
            print(f'error -> {genome}')
            continue

free_energy_dir = input('enter real sequence folding predictions: ')
free_energy_list = os.listdir(free_energy_dir)
try:
    free_energy_list.remove('.DS_Store')
except:
    pass

# cds and core at 37C
CDS = {}
for x in free_energy_list:
    name = x[:-11]

    with open(os.path.join(free_energy_dir, x)) as f:
        file = f.read().strip().split('>')
        CDS_ = []

        for line in file:
            gene = line.split('\n')

            try:
                data = gene[1]
            except:
                continue

            try:
                gene_CDS = string_parse(data)[2]
                if gene_CDS <= 0:
                    CDS_.append(gene_CDS)

            except:
                print(f'error at {x}')
                continue

        CDS_mean = np.mean(CDS_)
        CDS[name] = CDS_mean

core_dir = input('enter gene core folding dir: ').strip()
core_list = os.listdir(core_dir)
try:
    core_list.remove('.DS_Store')
except:
    pass

core = {}
for x in core_list:
    name = x[:-11]

    with open(os.path.join(core_dir, x)) as f:
        file = f.read().strip().split('>')
        CDS_ = []

        for line in file:
            gene = line.split('\n')

            try:
                data = gene[1]
            except:
                continue

            try:
                gene_CDS = string_parse(data)
                CDS_.append(gene_CDS)
            except:
                print(f'error at {x}')
                continue

        CDS_mean = np.mean(CDS_)
        core[name] = CDS_mean

with open(os.path.join(directory, f'metadata.csv'), 'w') as file:
    file.write(f"Genera,Growth_Temp,Genome_GC,CDS_room,Core_room\n")
    
    for genus in genera:
        try:
            file.write(f'{genus},{growth_dict[genus]},{gc_contents[genus]},{CDS[genus]},{core[genus]}\n')
        except:
            print(f'missing {genus}')
            continue
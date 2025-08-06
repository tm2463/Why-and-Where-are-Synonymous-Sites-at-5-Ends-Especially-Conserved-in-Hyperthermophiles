#list all core genes for building phylogeny
import os
import re
from Bio import SeqIO

input_dir = input('enter cds directory path: ').strip()
input_list = []

growth_file = input('enter growth file path: ').strip()
with open(growth_file) as f:
    filtered_list = f.read().strip().split('\n')

    for x in filtered_list:
        input_list.append(f'{x.split('\t')[0]}.txt')

input_list = [x for x in input_list]

#collect all protein names
protein_matrix = []
for x in input_list:
    with open(os.path.join(input_dir, x)) as f:
        protein_list = []

        try:
            for record in SeqIO.parse(f, 'fasta'):
                description = record.description.split('[')
                for header in description:
                    if header[:8] == 'protein=':
                        protein_list.append(header[8:-2])
            protein_matrix.append(protein_list)

        except Exception as e:
            print(x, e)

#create sets to remove duplicate protein names
genome_sets = [set(genus) for genus in protein_matrix]
#fine core genes across all genome sets
core_genes = set.intersection(*genome_sets)

#write output file
pwd = os.getcwd()
with open(os.path.join(pwd, 'core_genes.txt'), 'w') as output:
    for gene in core_genes:
        output.write(f'{gene}\n')

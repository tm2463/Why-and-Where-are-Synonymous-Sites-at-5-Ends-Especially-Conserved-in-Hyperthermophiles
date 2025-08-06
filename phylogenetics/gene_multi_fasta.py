#write multi fasta files for all genes identified in core_genes.py
import os
from Bio import SeqIO
from Bio.Seq import Seq

input_dir = input('enter cds directory path: ').strip()
gene_file = input('enter core gene file path: ').txt
growth_file = input('enter growth file path: ').txt
pwd = os.getcwd()
output_dir = os.path.join(pwd, 'multi_fasta')
os.makedirs(output_dir, exist_ok=True)

input_list = []
with open(growth_file) as f:
    filtered_list = f.read().strip().split('\n')
    for x in filtered_list:
        input_list.append(f'{x.split('\t')[0]}.txt')

try:
    input_list.remove('.DS_Store')
except:
    print('no ".DS_Store')

#read core genes file
with open(gene_file) as f:
    file_list = f.read().strip().split('\n')

#loop over all core genes
for gene in file_list:
    count = 0 #only allowing genes that feature once and only once in all species
    translations = []

    #write multi fasta file for each core gene from all species in the dataset
    try:
        for genus in input_list:
            with open(os.path.join(input_dir, genus)) as cds:
                for record in SeqIO.parse(cds, 'fasta'):

                    if f'[protein={gene}]' in record.description:
                        sequence = record.seq
                        translation = Seq(sequence).translate()
                        translations.append((genus[:-4], str(translation)))
                        count += 1 #quality control -> count printed during logging, should equal number of species
                        break

        try:
            #write file
            with open(os.path.join(output_dir, f'{gene}.txt'), 'w') as output:
                for genus, translation in translations:
                    output.write(f'>{genus}\n{translation}\n')

        except:
            #logging skipped genes
            print(f'skipping {gene} -> invalid file name')

    except Exception as e:
        #logging skipped genes
        print(f"skipping {gene} -> {e}")
        continue

    #logging
    print(f'{gene} -> {count}')
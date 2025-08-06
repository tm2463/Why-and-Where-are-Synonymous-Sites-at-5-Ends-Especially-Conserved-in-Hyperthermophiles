#concatenate MSA's and create partitions file for IQTREE
import os
from collections import defaultdict

input_dir = input('enter cleaned MSA directory path: ').strip()
input_list = os.listdir(input_dir)
try:
    input_list.remove('.DS_Store')
except:
    print('no ".DS_Store"')

concat = defaultdict(str)
partition_lines = []
current_start = 1

#loop over cleaned msa files
for n, i in enumerate(input_list):
    with open(os.path.join(input_dir, i)) as f:
        file = f.read().split('>')[1:]

        for x in file:
            entry = x.strip().split('\n')
            genus = entry[0]
            seqs = ''.join(entry[1:])

            if len(seqs) > 1000: #skip misaligned genes
                continue

            concat[genus] += seqs #add sequence string onto existing sequence in dictionary
            gene_length = len(seqs)

    #save lines for partition file
    current_end = current_start + gene_length - 1
    partition_lines.append(f'charset part{n+1} = {current_start}-{current_end};')
    current_start = current_end + 1

pwd = os.getcwd()

#write concatenated alignment file
with open(os.path.join(pwd, f'concat.fasta'), 'w') as output:
    for genus, seq in concat.items():
        output.write(f'>{genus}\n{seq}\n')

#write partitions file
with open(os.path.join(pwd, f'partitions.nex'), 'w') as p:
    p.write('#nexus\nbegin sets;\n') #nexus format

    for line in partition_lines:
        p.write(f'{line}\n')

    p.write(f'end;')

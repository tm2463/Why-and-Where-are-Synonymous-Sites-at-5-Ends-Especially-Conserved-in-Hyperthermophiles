import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
from tqdm import tqdm

# setup cds directory
input_dir = input('enter cds directory path: ').strip()
files = os.listdir(input_dir)
try:
    files.remove('.DS_Store')
except:
    print('no ".DS_Store"')
files = sorted(files)

# setup genome directory
genome_dir = input('enter genome directory path: ').strip()
genomes_list = os.listdir(genome_dir)
try:
    genomes_list.remove('.DS_Store')
except:
    print('no ".DS_Store')
genomes_list = sorted(genomes_list)

# setup output directories
os.makedirs('elongate_output', exist_ok=True)
pwd = os.getcwd()
output_dir = os.path.join(pwd, 'elongate_output')

elongate_dir = os.path.join(output_dir, 'elongate')
os.makedirs(elongate_dir, exist_ok=True)

intergenic_dir = os.path.join(output_dir, 'intergenic')
os.makedirs(intergenic_dir, exist_ok=True)

# extract coordinates from refseq cds files
ids, locations = [], []
for file in files:
    id, location = [], []

    with open(os.path.join(input_dir, file)) as f:
        for record in SeqIO.parse(f, 'fasta'):
            id.append(record.id)
            description = record.description
            match = re.search(r'\[location=(.*?)]', description)

            if match:
                location.append(match.group(1))

    ids.append(id)
    locations.append(location)

# load genome sequences
genomes = []
for genome in genomes_list:
    with open(os.path.join(genome_dir, genome)) as f:
        lines = f.read().split('\n')
        sequence = ''.join(lines[1:])
        genomes.append(sequence)

# write output files (5' seqs (cds -30nt to +60nt relative to start & intergenome)
for genome, file, id, location in tqdm(zip(genomes, files, ids, locations), desc='processing'):
    starts, ends = [], []

    with open(os.path.join(elongate_dir, file), 'w') as e_out:
        for seq_id_entry, coord in zip(id, location):
            complement = coord.startswith('complement')  # complement requires different processing

            if complement:
                coords = coord[11:-1]
                coords = coords.split('..')
            else:
                coords = coord.split('..')

            try:  # if no coordinates are found, skip
                start = int(coords[0])
                end = int(coords[1])
            except:
                continue

            starts.append(start)
            ends.append(end)

            if len(starts) < 2:
                continue

            if start > starts[-2]:  # plasmid skipping logic
                if complement:
                    end += 30
                    fragment = str(Seq(genome[start + 1:end]).reverse_complement())
                else:
                    start = max(0, start - 31)
                    fragment = genome[start:end]

                e_out.write(f'>{seq_id_entry}\n{fragment[:90]}\n')
            else:
                break  # once chromosome has been covered, skip to next sample

    # write intergenome to file
    with open(os.path.join(intergenic_dir, file), 'w') as i_out:
        ends.insert(0, 0)
        ends = ends[:-1]
        intergenic = []

        for start, end in zip(starts, ends):
            inter_sec = genome[end:start]
            intergenic.append(inter_sec)

        intergenic = ''.join([x for x in intergenic])
        i_out.write(f'{intergenic}')

import numpy as np
import pandas as pd
import os
import requests
import gzip
import io

os.makedirs('results', exist_ok=True)
os.chdir('results')
path = os.getcwd()
summary = os.path.join(path, 'assembly_summary.txt')

#download assembly summary file from NCBI FTP
summary_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'

ncbi = requests.get(summary_url)

with open(summary, 'wb') as f:
	f.write(ncbi.content)

df = pd.read_csv(summary, sep='\t', skiprows=1)
data = df[['#assembly_accession', 'taxid', 'organism_name', 'assembly_level', 'ftp_path', 'group', 'genome_size', 'genome_size_ungapped', 'gc_percent', 'total_gene_count', 'protein_coding_gene_count', 'non_coding_gene_count']]

"""
filter conditions:
assembly level = complete genome or chromosome
1 genome per genus as species are not independent data points (assembly with least gaps in case of extra assemblies)
"""

genus = np.array(data['organism_name'].str.split().str[0])
data['genus'] = genus
filtered_data = data[data['assembly_level'] == 'Complete Genome']
filtered_data = filtered_data.groupby('genus', as_index=False).first()
filtered_data.to_csv('filtered_assemblies.csv', index=False)

#make directories and store file suffix
os.makedirs('cds', exist_ok=True)
cds_dir = os.path.join(path, 'cds')
cds_file = '_cds_from_genomic.fna.gz'

os.makedirs('genome', exist_ok=True)
genome_dir = os.path.join(path, 'genome')
genome_file = '_genomic.fna.gz'

data = pd.read_csv('filtered_assemblies.csv')
url_dict = dict(zip(data['genus'], data['ftp_path']))

#download data
for n, genus in enumerate(url_dict, 1):
	try:
		#create urls leading to data
		link = url_dict[genus]
		accession = link.split('/')[-1]

		cds_url = f"{link}/{accession}{cds_file}"
		genome_url = f"{link}/{accession}{genome_file}"

		#get data
		cds_data = requests.get(cds_url, timeout=120)
		genome_data = requests.get(genome_url, timeout=120)

		#write cds file
		with gzip.open(io.BytesIO(cds_data.content), 'rb') as input:
			with open(os.path.join(cds_dir, f'{genus}.txt'), 'wb') as output:
				output.write(input.read())

		#write genome file
		with gzip.open(io.BytesIO(genome_data.content), 'rb') as input:
			with open(os.path.join(genome_dir, f'{genus}.txt'), 'wb') as output:
				output.write(input.read())

		print(f'{genus} ({n}/{len(url_dict)}) -> complete')
		
	except:
		print(f'error downloading {genus}')

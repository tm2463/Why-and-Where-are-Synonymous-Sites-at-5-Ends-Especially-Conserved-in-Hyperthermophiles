#get references for salmon
import pandas as pd

assembly_sum = input('enter assembly summary file path: ').strip()
master = input('enter master.tsv file path: ').strip()

with open(assembly_sum) as f:
    assem_sum = pd.read_csv(f)
    ref_df = assem_sum[['organism_name', 'refseq_category', 'ftp_path']]

with open(master) as f:
    names = pd.read_csv(f, sep='\t', header=None)
    names = names[1].tolist()

names = [name.replace('_', ' ') for name in names]

ref_df = ref_df[ref_df['organism_name'].isin(names)]
ref_df = ref_df[ref_df['refseq_category'] == 'reference genome']

ref_df.to_csv('reference.csv', index=False)

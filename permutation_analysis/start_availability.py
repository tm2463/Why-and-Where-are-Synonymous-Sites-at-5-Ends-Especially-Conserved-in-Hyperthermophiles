#modified version of real_fold.py
import RNA
import os
from concurrent.futures import ProcessPoolExecutor, as_completed


#fold all genes in the genus file
def fold_genus_file(genus_file, fold_temp, growth_dict, path, results_dir):
    genus_name = genus_file[:-4]
	
	#control whether folding at OGT or set temp
    if fold_temp == 'ogt'.lower():
        temp = growth_dict[genus_name]
    else:
        temp = int(fold_temp)

    RNA.cvar.temperature = temp #set folding temp, if OGT, temp will vary with each genome

	#load genes
    with open(os.path.join(path, genus_file)) as f:
        genes = f.read().split('>')[1:]
	
	#define out path
    out_path = os.path.join(results_dir, f'{genus_name}_energy.txt')
    with open(out_path, 'w') as out:
        for id, gene_entry in enumerate(genes):
            lines = gene_entry.strip().split('\n')
            seq = lines[1].strip()
            fold = RNA.fold(seq[15:45]) #only fold start of seq

            if float(fold[1]) > 1: #avoid strange behaviour at end of files
                print(f'{genus_name} skipping gene: {id}')
                continue
                
			#write the structure prediction, only start codon is written to file
            out.write(f'>{id}_{genus_file}\n{fold[0][:3]}\n')

    print(f'finished {genus_name}')

if __name__ == '__main__':
    path = input('enter elongated sequences directory path: ').strip()
    growth_file = input('enter growth temps file path: ').strip()
    fold_temp = input('temp temp selection: if OGT, enter "ogt", otherwise enter temp value (°C): ')
    
    pwd = os.getcwd()
    output_dir = os.path.join(pwd, f'real_fold_{fold_temp}')
    results_dir = os.path.join(output_dir, f'results_{str(fold_temp)}')
    os.makedirs(results_dir, exist_ok=True)

    # load growth temps file
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

    real_seq = os.listdir(path)
    if '.DS_Store' in real_seq:
        real_seq.remove('.DS_Store')

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        for genus in real_seq:
            executor.submit(fold_genus_file, genus, fold_temp, growth_dict, path, results_dir)

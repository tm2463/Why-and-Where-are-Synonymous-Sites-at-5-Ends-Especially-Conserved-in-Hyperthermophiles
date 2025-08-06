#output directories/files can be renamed for sake of organisation
import RNA
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor, as_completed


def average_permutation(genus_path, temp):
    genes = os.listdir(genus_path)
    gene_folds = []

    RNA.cvar.temperature = temp

    for x in genes:
        gene_mean = []

        with open(os.path.join(genus_path, x)) as f:
            file = f.read().split('>')[1:]
            fold_1, fold_2, fold_3, fold_4, fold_5 = [], [], [], [], []

            for permutation in file:
                seq = permutation.split('\n')[1]
                fold_1.append(RNA.fold(seq[:30])[1])
                fold_2.append(RNA.fold(seq[15:45])[1])
                fold_3.append(RNA.fold(seq[30:60])[1])
                fold_4.append(RNA.fold(seq[45:75])[1])
                fold_5.append(RNA.fold(seq[60:90])[1])

        gene_mean.append(float(np.mean(fold_1)))
        gene_mean.append(float(np.mean(fold_2)))
        gene_mean.append(float(np.mean(fold_3)))
        gene_mean.append(float(np.mean(fold_4)))
        gene_mean.append(float(np.mean(fold_5)))

        gene_folds.append(gene_mean)
    return gene_folds


def process_genus(genus, input_dir, output_dir, growth_temp_file):
    print(f'starting processing {genus}')

    with open(growth_temp_file) as f:
        growth_temps = f.read().strip().split('\n')

    growth_dict = {}
    for x in growth_temps:
        try:
            y = x.split('\t')
            g = y[0]
            temp = y[1]
            growth_dict[g] = float(temp)
        except:
            print(f'missing temps for {g}')
            return

    temp = growth_dict[genus]
    genus_path = os.path.join(input_dir, genus)
    out_file_path = os.path.join(output_dir, f'{genus}_permutation.txt')

    try:
        with open(out_file_path, 'w') as file:
            for x, gene_folds in enumerate(average_permutation(genus_path, temp)):
                file.write(f'>{x}_{genus}\n')
                file.write(f'{gene_folds}\n')
    except:
        print(f'error folding {genus}')


def main():
    input_dir = input('enter path to permutations directory: ').strip()
    os.makedirs('mean_permutations_output')
    pwd = os.getcwd()
    output_dir = os.path.join(pwd, 'mean_permutations_output')
    growth_temp_file = input('enter growth temp file path: ').strip()

    permutation_results = os.path.join(output_dir, 'permutation_results')
    os.makedirs(permutation_results, exist_ok=True)

    genus_list = [f for f in os.listdir(input_dir) and f != '.DS_Store']

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        for genus in genus_list:
            executor.submit(process_genus, genus, input_dir, permutation_results, growth_temp_file)


if __name__ == '__main__':
    main()

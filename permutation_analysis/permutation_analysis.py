import os
import numpy as np

real_results = input('enter real folding predictions directory path: ').strip()
permutations_results = input('enter permuted folding predictions directory path: ').strip()
growth_file = input('enter growth file path: ').strip()
genome_dir = input('enter synonymous shuffle sequence output directory path: ').strip()

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

genome_list = os.listdir(genome_dir)
if '.DS_Store' in genome_list:
    genome_list.remove('.DS_Store')
    
mesophiles, thermophiles, hyperthermophiles = [], [], []
for x in genome_list:
    try:
        temp = growth_dict[x]
        if temp < 50:
            mesophiles.append(x)
        elif 50 <= temp < 80:
            thermophiles.append(x)
        elif temp > 80:
            hyperthermophiles.append(x)
    except:
        continue


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


def get_z_sem(group):
    real_files = os.listdir(real_results)
    perm_files = os.listdir(permutations_results)

    real = sorted([f.replace('_energy.txt', '') for f in real_files])
    perm = sorted([f.replace('_permutation.txt', '') for f in perm_files])

    genus_list = [x for x in real if x in perm]
    genus_list = [x for x in genus_list if x in group]

    filtered_real = [f'{x}_energy.txt' for x in genus_list]
    filtered_perm = [f'{x}_permutation.txt' for x in genus_list]

    real_lists = []
    perm_lists = []
    for x in range(len(filtered_real)):
        try:
            real_file = filtered_real[x]
            perm_file = filtered_perm[x]

            with open(os.path.join(real_results, real_file)) as real_f:
                real_content = real_f.read().split('>')[1:]
                real_fold_data = [string_parse(entry.split('\n')[1]) for entry in real_content]

            with open(os.path.join(permutations_results, perm_file)) as perm_f:
                perm_content = perm_f.read().split('>')[1:]
                perm_fold_data = [string_parse(entry.split('\n')[1]) for entry in perm_content]

            if real_fold_data and perm_fold_data:
                real_lists.append(real_fold_data)
                perm_lists.append(perm_fold_data)

        except:
            continue

    real_means = np.zeros((len(real_lists), 5), dtype=float)
    perm_means = np.zeros((len(perm_lists), 5), dtype=float)
    
    for n, (x, y) in enumerate(zip(real_lists, perm_lists)):
        real_arr = np.array(x)
        perm_arr = np.array(y)

        real_means[n] = np.mean(real_arr, axis=0)
        perm_means[n] = np.mean(perm_arr, axis=0)

    real_mean = np.mean(real_means, axis=0)
    perm_mean = np.mean(perm_means, axis=0)
    perm_stdev = np.std(perm_means, axis=0)
    group_size = len(real_lists)

    z_score = [(rm - pm)/ps for rm, pm, ps in zip(real_mean, perm_mean, perm_stdev)]
    sem = [ps/np.sqrt(group_size) for ps in perm_stdev]

    return z_score, sem


meso_z, meso_sem = get_z_sem(mesophiles)
thermo_z, thermo_sem = get_z_sem(thermophiles)
hyper_z, hyper_sem = get_z_sem(hyperthermophiles)

pwd = os.getcwd()
os.makedirs(os.path.join(pwd, 'permutations_ogt'), exist_ok=True)
out_dir = os.path.join(pwd, 'permutations_ogt')

meso = np.column_stack((meso_z, meso_sem))
np.savetxt(os.path.join(out_dir, 'meso.csv'), meso, delimiter=',')

thermo = np.column_stack((thermo_z, thermo_sem))
np.savetxt(os.path.join(out_dir, 'thermo.csv'), thermo, delimiter=',')

hyper = np.column_stack((hyper_z, hyper_sem))
np.savetxt(os.path.join(out_dir, 'hyper.csv'), hyper, delimiter=',')

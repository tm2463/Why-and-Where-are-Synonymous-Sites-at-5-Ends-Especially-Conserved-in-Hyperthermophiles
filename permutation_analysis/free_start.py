import os
import numpy as np

growth_file = input('enter growth file path: ').strip()
real_results = input('enter real folding predictions directory path: ').strip()
permutations_results = input('enter permuted folding predictions directory path: ').strip()
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

real_list = os.listdir(real_results)
try:
    real_list.remove('.DS_Store')
except:
    pass

permuted_list = os.listdir(permutations_results)
try:
    permuted_list.remove('.DS_Store')
except:
    pass

real_start_dict = {}
for x in real_list:
    with open(os.path.join(real_results, f'{x}')) as f:
        file = f.read().strip().split('>')
        starts = []

        for gene in file:
            try:
                seq = gene.split('\n')[1]
                starts.append(seq)
            except:
                continue

        starts = ''.join([s for s in starts])
        length = len(starts)
        empty_count = starts.count('.')
        free_codon = empty_count/length
        real_start_dict[x[:-11]] = free_codon

perm_start_dict = {}
for x in permuted_list:
    with open(os.path.join(permutations_results, f'{x}')) as f:
        file = f.read().strip().split('>')
        starts = []

        for gene in file:
            try:
                seq = gene.split('\n')[1]
                starts.append(float(seq))
            except:
                continue

        mean_starts = np.mean(starts)
        perm_start_dict[x[:-16]] = mean_starts

pwd = os.getcwd()
t_dir = os.path.join(pwd, 'ttest')
os.makedirs(t_dir, exist_ok=True)

meso_real = np.array([real_start_dict[x] for x in mesophiles if x in real_start_dict])
meso_perm = np.array([perm_start_dict[x] for x in mesophiles if x in real_start_dict])
thermo_real = np.array([real_start_dict[x] for x in thermophiles if x in real_start_dict])
thermo_perm = np.array([perm_start_dict[x] for x in thermophiles if x in real_start_dict])
hyper_real = np.array([real_start_dict[x] for x in hyperthermophiles if x in real_start_dict])
hyper_perm = np.array([perm_start_dict[x] for x in hyperthermophiles if x in real_start_dict])

meso = np.column_stack((meso_real, meso_perm))
thermo = np.column_stack((thermo_real, thermo_perm))
hyper = np.column_stack((hyper_real, hyper_perm))

np.savetxt(os.path.join(t_dir, 'meso.csv'), meso, delimiter=',')
np.savetxt(os.path.join(t_dir, 'thermo.csv'), thermo, delimiter=',')
np.savetxt(os.path.join(t_dir, 'hyper.csv'), hyper, delimiter=',')

del_dir = os.path.join(pwd, 'delta_start')
os.makedirs(del_dir, exist_ok=True)

delta_meso = np.array([real - perm for real, perm in zip(meso_real, meso_perm)])
delta_thermo = np.array([real - perm for real, perm in zip(thermo_real, thermo_perm)])
delta_hyper = np.array([real - perm for real, perm in zip(hyper_real, hyper_perm)])

np.savetxt(os.path.join(del_dir, 'meso.csv'), delta_meso, delimiter=',')
np.savetxt(os.path.join(del_dir, 'thermo.csv'), delta_thermo, delimiter=',')
np.savetxt(os.path.join(del_dir, 'hyper.csv'), delta_hyper, delimiter=',')

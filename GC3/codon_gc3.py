import numpy as np
import os

growth_file = input('enter growth file path: ').strip()
genome_dir = input('enter genome directory path: ').strip()
cds_dir = input('enter cds directory path: ').strip()

#load growth temps file
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
try:
    genome_list.remove('.DS_Store')
except:
    print('no ".DS_Store"')

mesophiles, thermophiles, hyperthermophiles = [], [], []
for x in genome_list:
    try:
        if growth_dict[x[:-4]] < 50:
            mesophiles.append(x)
        elif growth_dict[x[:-4]] > 50 and growth_dict[x[:-4]] < 80:
            thermophiles.append(x)
        elif growth_dict[x[:-4]] > 80:
            hyperthermophiles.append(x)
    except:
        continue

genome_list = []
genome_list.append(mesophiles)
genome_list.append(thermophiles)
genome_list.append(hyperthermophiles)

cds_list = os.listdir(cds_dir)
try:
    cds_list.remove('.DS_Store')
except:
    print('no ".DS_Store"')


mesophile_cds, thermophile_cds, hyperthermophile_cds = [], [], []
for group, new_list in zip(genome_list, [mesophile_cds, thermophile_cds, hyperthermophile_cds]):
    for genus in group:
        with open(os.path.join(cds_dir, f'{genus}')) as f:
            data = f.read().split('>')
            cds = []

            for x in data:
                seq = ''.join(x.split('\n')[1:])

                if seq:
                    cds.append(seq[:45])

        new_list.append(cds)


def gc3_slicer(group):
    nt_dict = {'C': 1, 'G': 1, 'A': 0, 'T':0}
    group_arr = []

    for genera in group:
        string_list = []

        for gene in genera:
            third_pos = [nt_dict[gene[x+2]] for x in range(0, len(gene), 3)]
            third_pos = third_pos[1:]
            string_list.append(third_pos)

        string_arr = np.array(string_list)
        genera_arr = np.nanmean(string_arr, axis=0)
        group_arr.append(genera_arr)

    whole_arr_mean = np.nanmean(np.array(group_arr))
    position_means = np.nanmean(group_arr, axis=0)
    position_stdev = np.std(group_arr, axis=0)

    return whole_arr_mean, position_means, position_stdev

meso_arr_mean, meso_means, meso_stdevs = gc3_slicer(mesophile_cds)
thermo_arr_mean, thermo_means, thermo_stdevs = gc3_slicer(thermophile_cds)
hyper_arr_mean, hyper_means, hyper_stdevs = gc3_slicer(hyperthermophile_cds)

#calculate z scores and sem
meso_size = len(mesophile_cds)
thermo_size = len(thermophile_cds)
hyper_size = len(hyperthermophile_cds)

meso_z = (meso_means - meso_arr_mean) / meso_stdevs
meso_sem = meso_stdevs / np.sqrt(meso_size)

thermo_z = (thermo_means - thermo_arr_mean) / thermo_stdevs
thermo_sem = thermo_stdevs / np.sqrt(thermo_size)

hyper_z = (hyper_means - hyper_arr_mean) / hyper_stdevs
hyper_sem = thermo_stdevs / np.sqrt(hyper_size)

pwd = os.getcwd()
os.makedirs(os.path.join(pwd, 'gc3'), exist_ok=True)
out_dir = os.path.join(pwd, 'gc3')

meso = np.column_stack((meso_z, meso_sem))
np.savetxt(os.path.join(out_dir, 'meso.csv'), meso, delimiter=',')

thermo = np.column_stack((thermo_z, thermo_sem))
np.savetxt(os.path.join(out_dir, 'thermo.csv'), thermo, delimiter=',')

hyper = np.column_stack((hyper_z, hyper_sem))
np.savetxt(os.path.join(out_dir, 'hyper.csv'), hyper, delimiter=',')
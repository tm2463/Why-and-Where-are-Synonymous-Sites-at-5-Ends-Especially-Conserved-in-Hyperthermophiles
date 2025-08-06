import os
import numpy as np


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


directory = os.getcwd()
growth_file = input('enter predicted OGT file path: ').strip()

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

#ensure real folds for pgls are saved in folders results_40, results_50, etc...
pgls_temp = input('enter pgls real fold directory path: ').strip()

forty = os.path.join(pgls_temp, 'results_40')
forty_list = [f for f in os.listdir(forty) if f != '.DS_Store']
fifty = os.path.join(pgls_temp, 'results_50')
fifty_list = [f for f in os.listdir(fifty) if f != '.DS_Store']
sixty = os.path.join(pgls_temp, 'results_60')
sixty_list = [f for f in os.listdir(sixty) if f != '.DS_Store']
seventy = os.path.join(pgls_temp, 'results_70')
seventy_list = [f for f in os.listdir(seventy) if f != '.DS_Store']
eighty = os.path.join(pgls_temp, 'results_80')
eighty_list = [f for f in os.listdir(eighty) if f != '.DS_Store']
ninety = os.path.join(pgls_temp, 'results_90')
ninety_list = [f for f in os.listdir(ninety) if f != '.DS_Store']


def extract_means(file_list, file_path):
    func_dict = {}

    for genera in file_list:
        with open(os.path.join(file_path, genera)) as f:
            file = f.read().strip().split('>')
            values = []

            for gene in file:
                try:
                    data = gene.split('\n')[1]
                except:
                    continue

                values.append(float(data))

            genus_mean = np.mean(values)
            func_dict[genera[:-11]] = genus_mean

    return func_dict


prime_forty = extract_means(forty_list, forty)
prime_fifty = extract_means(fifty_list, fifty)
prime_sixty = extract_means(sixty_list, sixty)
prime_seventy = extract_means(seventy_list, seventy)
prime_eighty = extract_means(eighty_list, eighty)
prime_ninety = extract_means(ninety_list, ninety)

genome_dir = input('enter genome directory path: ').strip()
genome_list = os.listdir(genome_dir)
if '.DS_Store' in genome_list:
    genome_list.remove('.DS_Store')

genera = [x[:-4] for x in genome_list]

with open(os.path.join(directory, f'metadata2.csv'), 'w') as f:
    f.write(f"Genera,Growth_Temp,40,50,60,70,80,90\n")

    for genus in genera:
        try:
            f.write(f'{genus},{growth_dict[genus]},{prime_forty[genus]},{prime_fifty[genus]},{prime_sixty[genus]},{prime_seventy[genus]},{prime_eighty[genus]},{prime_ninety[genus]}\n')
        except:
            print(f'missing {genus}')
            continue
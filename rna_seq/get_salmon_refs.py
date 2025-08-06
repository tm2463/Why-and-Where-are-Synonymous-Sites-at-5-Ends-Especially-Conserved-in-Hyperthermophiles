#downloads reference and write file in correct folder for mapping with salmon
import os
import pandas as pd
import requests
import gzip
import io

ref = input('enter reference.csv file path: ').strip()

with open(ref) as f:
    df = pd.read_csv(f, header=None, skiprows=1)

genera = df[0].tolist()
links = df[2].tolist()

type = '_cds_from_genomic.fna.gz'
for x, url in zip(genera, links):

    try:
        os.chdir(os.path.join(ref, 'rnaseq_data', x))
    except:
        print(f'no directory -> {x}')
        continue

    filename = '/' + url.split('/')[-1] + type
    link = url + filename

    try:
        response = requests.get(link, timeout=120)
        with gzip.open(io.BytesIO(response.content), 'rb') as input:
            with open(f'{filename[:-3]}', 'wb') as output:
                output.write(input.read())
                print(f'downloaded {x} reference')

    except:
        print(f'timeout at {x}')

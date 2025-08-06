# takes complete genome sequences scraped from web_scraper.py and returns predicted optimal growth temp
# OGT predicted from http://www.orgene.net/CnnPOGTP/predictor2.php
import requests
import os
from bs4 import BeautifulSoup

genome_dir = input('enter genome directory path: ').strip()
output_dir = input('enter output file path: ').strip()

url = "http://www.orgene.net/CnnPOGTP/predictor2.php"
file_names = [x for x in os.listdir(genome_dir)]

for file in file_names:
    file_path = os.path.join(genome_dir, file)

    with open(file_path, 'rb') as f:
        try:
            files = {'genome_fasta_file': f}
            print(f'sending {file} to server')
            response = requests.post(url, files=files, timeout=300)

            if response.ok:
                print(f'{file} response ok')
                soup = BeautifulSoup(response.text, 'html.parser')
                temp = soup.find('font', class_='Introduction_font3')
                text = temp.text.strip().replace('&nbsp', '').replace('℃', '')
                entry = f'{file[:-4]}\t{text}'
            else:
                entry = f'{file} invalid\n'
        except:
            print(f'{file} -> failed')

    with open(os.path.join(output_dir, 'growth_temps.txt'), 'a') as out:
        out.write(entry)
        print(entry)
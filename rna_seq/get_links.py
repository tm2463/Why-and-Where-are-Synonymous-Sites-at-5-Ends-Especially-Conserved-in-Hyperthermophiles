import pandas as pd
import os
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

path = os.getcwd()
summary = os.path.join(path, 'assembly_summary.txt')

summary_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'

ncbi = requests.get(summary_url)

with open(summary, 'wb') as f:
    f.write(ncbi.content)

summary_df = pd.read_csv(summary, sep='\t', skiprows=1)
summary_df = summary_df[['species_taxid', 'organism_name', 'ftp_path']]
summary_df = summary_df.drop_duplicates(subset='species_taxid', keep='first')
species_id = summary_df['species_taxid'].tolist()
searches = [f'https://www.ncbi.nlm.nih.gov/sra/?term=txid{x}%5BOrganism%3Aexp%5D+RNA-Seq' for x in species_id]

driver = webdriver.Chrome()
timeout = 2.5

links = []
for i, url in enumerate(searches):
    print(f"{i+1}/{len(searches)}: {url}")

    try:
        driver.get(url)
        wait = WebDriverWait(driver, timeout=timeout)
        result = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, '.rprt .title a')))
        link = result.get_attribute('href')
        print(f"valid result -> {link}")
        links.append(link)

    except:
        print(f'no valid results')

driver.quit()

with open(os.path.join(path, 'links.txt')) as output:
    for link in links:
        output.write(f'{link}\n')

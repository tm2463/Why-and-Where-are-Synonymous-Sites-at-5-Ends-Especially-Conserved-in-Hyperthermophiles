from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import requests

sra = input('enter path to links file: ').strip()

with open(sra) as f:
    links = f.read().strip().split('\n')

driver = webdriver.Chrome()
timeout = 10

runs, species = [], []
for i, url in enumerate(links):
    print(f"{i + 1}/{len(links)}: {url}")

    accession = url.strip('/').split('/')[-1]

    driver.get(url)
    wait = WebDriverWait(driver, timeout=timeout)
    result = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, 'tr a')))
    link = result.get_attribute('href')
    run = link.split('=')[-1]
    runs.append(run)

    html = requests.get(url)
    soup = BeautifulSoup(html.text, 'html.parser')
    items = soup.find_all('div', class_='sra-full-data')

    for item in items:
        if 'Organism' in item.text:
            organism = item.text.split(': ')[-1]
            organism = organism.replace(' ', '_')

    species.append(organism)

with open('master.tsv', 'w') as f:
    for run, spec in zip(runs, species):
        f.write(f"{run}\t{spec}\n")

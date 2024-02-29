import time
import os
import json

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC


def upload(img_path):
    img_path = "C:\\pycharmpro\\chemModule\\structuralFormuls\\US20230136730A1-20230504-C00209.png"
    select_xpath = 'body > center:nth-child(1) > form:nth-child(3) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > input:nth-child(3)'
    submit_xpath = '#b_upload'
    getSmiles = '#b_getsmiles'
    #clear_xpath = '/html/body/center/form/table/tbody/tr[3]/td[1]/center/input[2]'

    #firefox.find_element(By.XPATH(clear_xpath)).click()
    firefox.find_element(By.CSS_SELECTOR, select_xpath).send_keys(img_path)
    firefox.find_element(By.CSS_SELECTOR, submit_xpath).click()
    #time.sleep(20)
    #firefox.find_element(By.CSS_SELECTOR, getSmiles).click()


def get_information():
    get_smiles_xpath = '#b_getsmiles'
    smiles_xpath = 'body > center:nth-child(1) > form:nth-child(3) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(2) > input:nth-child(1)'

    firefox.find_element(By.CSS_SELECTOR, get_smiles_xpath).click()
    text = firefox.find_element(By.CSS_SELECTOR, smiles_xpath).get_attribute("value")

    return text


def main():
    global firefox
    #chrome_driver_path = os.path.abspath("chromedriver.exe")  # Получаем абсолютный путь к chromedriver.exe

    options1 = Options()
    options1.add_argument('--headless=new')
    firefox = webdriver.Chrome(options1)
    firefox.get('https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi')
    wait = WebDriverWait(firefox, 20)

    img_folder = 'structuralFormuls'
    imgs76 = os.listdir(img_folder)
    smiles_list = {}
    for img in imgs76:
        upload(img_folder.replace('/', '\\') + '\\' + img)
        time.sleep(7)
        try:
            tmp_text = get_information()
            firefox.save_screenshot('res/' + img.rstrip('.jpg') + '.png')
            smiles_list[img] = tmp_text
        except:
            smiles_list[img] = 'Sorry, no structures found'

        print(smiles_list[img])
    firefox.quit()

    with open('result.json', 'w') as fp:
        json.dump(smiles_list, fp)


if __name__ == '__main__':
    main()

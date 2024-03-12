import time
import os
import json

from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By

def upload(img_path):
    img_path_root = os.path.join(os.getcwd())
    img_path = img_path_root + os.sep + img_path
    select_xpath = 'body > center:nth-child(1) > form:nth-child(3) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > input:nth-child(3)'
    submit_xpath = '#b_upload'

    driver.find_element(By.CSS_SELECTOR, select_xpath).send_keys(img_path)
    driver.find_element(By.CSS_SELECTOR, submit_xpath).click()

def get_information():
    get_smiles_xpath = '#b_getsmiles'
    smiles_xpath = 'body > center:nth-child(1) > form:nth-child(3) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(2) > input:nth-child(1)'

    driver.find_element(By.CSS_SELECTOR, get_smiles_xpath).click()
    text = driver.find_element(By.CSS_SELECTOR, smiles_xpath).get_attribute("value")

    return text

def IMG2SMILES(img_path1):
    # основная процедура
    global driver
    chrome_options = webdriver.ChromeOptions()

    # Включаем использование DevTools
    chrome_options.add_experimental_option("excludeSwitches", ["enable-automation"])
    chrome_options.add_experimental_option("useAutomationExtension", False)
    chrome_options.add_argument("--disable-blink-features=AutomationControlled")
    chrome_options.add_argument("--disable-infobars")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--disable-browser-side-navigation")
    chrome_options.add_argument("--disable-gpu")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-notifications")
    chrome_options.add_argument("--headless")

    # Запускаем Chrome с DevTools Protocol (CDP)
    driver = webdriver.Edge(options=chrome_options)
    driver.get('https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi')
    wait = WebDriverWait(driver, 0)

    img_folder = 'structuralFormuls'
    imgs76 = os.listdir(img_folder)
    smiles_list = {}


    # Вызываем функцию upload для загрузки каждого изображения
    upload(img_path1)
    time.sleep(4)
    try:
        tmp_text = get_information()
        smiles_list[img_path1] = tmp_text
    except:
        smiles_list[img_path1] = 'Sorry, no structures found'
    print(smiles_list[img_path1])



    driver.quit()

    with open('result.json', 'w') as fp:
        json.dump(smiles_list, fp)

    return smiles_list[img_path1]

IMG2SMILES("Patent_images\\US7754717B2\\US7754717B2_image_32.png")
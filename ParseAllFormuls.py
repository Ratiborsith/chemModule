import sqlite3
import time
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
import urllib.request
import os
from selenium.common.exceptions import TimeoutException
import requests
import re
from bs4 import BeautifulSoup

def get_patent_id(cursor, patent_number):
    cursor.execute("SELECT id FROM patents WHERE PatentNumber = ?", (patent_number,))
    result = cursor.fetchone()
    if result:
        return result[0]
    else:
        return None


def download_image(url, filename, max_retries=3):
    attempt = 0
    while attempt < max_retries:
        try:
            urllib.request.urlretrieve(url, filename)
            print(f"Image saved: {filename}")
            return
        except urllib.error.URLError as e:
            print(f"URL error occurred when trying to retrieve {url}: {e}")
            attempt += 1
            time.sleep(2)  # Wait 2 seconds before retrying
        except Exception as e:
            print(f"An error occurred when trying to retrieve {url}: {e}")
            attempt += 1
            time.sleep(2)
    print(f"Failed to download image after {max_retries} attempts.")


# Function to fetch webpage content
def fetch_webpage(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print("Failed to fetch webpage")
        return None

# Function to extract chemical formulas from HTML content
def extract_chemical_info(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    chemical_info = {}
    # Find all li elements with itemprop="match" attribute
    li_elements = soup.find_all('li', itemprop='match')
    for li in li_elements:
        # Find span elements with relevant itemprop values
        domain_span = li.find('span', itemprop='domain')
        name_span = li.find('span', itemprop='name')
        smiles_span = li.find('span', itemprop='smiles')
        # Check if all relevant information is present and domain is "Chemical compound"
        if domain_span and domain_span.text.strip() == "Chemical compound" and name_span and smiles_span:
            domain = domain_span.text.strip()
            name = name_span.text.strip()
            smiles = smiles_span.text.strip()
            # Use regular expressions to identify types of formulas
            if re.match(r"^[^=;]*$", name):  # UIPAC
                chemical_info.setdefault("UIPAC", []).append(name)
            elif re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#]+$', name):  # SMILES
                chemical_info.setdefault("SMILES", []).append(smiles)
            else:  # compoundName
                chemical_info.setdefault("compoundName", []).append(name)
    return chemical_info

# Main function
def parse_patent(patent_number_url):

    # сначала парсим все uipac - представления, а также текстовые, если встретятся
    url = f"https://patents.google.com/patent/{patent_number_url}"

    uipac_list = []
    webpage_content = fetch_webpage(url)
    if webpage_content:
        chemical_info = extract_chemical_info(webpage_content)
        # Print by the specified order
        i = 0
        for type_name in ["SMILES", "UIPAC", "compoundName"]:

            if type_name == "UIPAC":
                print(type_name + ":")
                for formula in chemical_info[type_name]:
                    print(formula)
                    uipac_list.append(formula)
                    i = i + 1


        print(f"Количество найденных формул: {i}")


    # Далее парсим все структурные представления формул, которые есть в патенте

    url = "https://patents.google.com/"
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument("--headless")  # Runs Chrome in headless mode.
    chrome_options.add_argument(
        "--disable-gpu")  # Disables GPU hardware acceleration. If software renderer is not in place, then the headless browser will not launch on Windows.
    chrome_options.add_argument(
        "--no-sandbox")  # Bypass OS security model. This can be required if running as root on Linux.
    chrome_options.add_argument("--disable-dev-shm-usage")  # Overcome limited resource problems.

    service = Service()
    browser = webdriver.Chrome(service=service, options=chrome_options)
    browser.get(url)

    patent_name = patent_number_url

    main_directory = "Patent_images"
    if not os.path.exists(main_directory):
        os.makedirs(main_directory)



    # Assuming 'patent_numbers' is a list of patent numbers and 'url' is the URL to start with
    #for patent in patent_numbers:
    browser.get(url)
    # Wait for the search box to be present and interact with it
    search_box = WebDriverWait(browser, 10).until(
        EC.presence_of_element_located((By.ID, "searchInput")))
    search_box.clear()
    search_box.send_keys(patent_name)   # здесь осуществляется поиск по имени
    search_button = WebDriverWait(browser, 10).until(
        EC.element_to_be_clickable((By.ID, "searchButton")))
    search_button.click()

    # Wait for the image carousel to be present
    try:
        image_carousel = WebDriverWait(browser, 10).until(EC.presence_of_element_located((By.ID, 'figures')))
        images = image_carousel.find_elements(By.TAG_NAME, 'img')
        flag = 0
    except TimeoutException:
        try:
            images = WebDriverWait(browser, 10).until(
                EC.presence_of_all_elements_located((By.CSS_SELECTOR, 'img[alt^="Figure"]')))
            flag = 1
        except TimeoutException:
            # No images found after all attempts
            print(f"В патенте {patent_name} не найдено структурных представлений.")


    print(f"Всего для патента {patent_name} найдено изображений: {len(images)}")
    current_patent_directory = os.path.join(main_directory, patent_name)

    if not os.path.exists(current_patent_directory):
        os.makedirs(current_patent_directory)
    structural_formul_list = []
    for index, img in enumerate(images, start=1):
        try:
            if flag == 1:
                image_src = img.get_attribute('src')
                image_filename = os.path.join(current_patent_directory, f"{patent_name}_image_{index}.png")
                download_image(image_src, image_filename)
                # print(f"Image saved: {image_filename}")
                structural_formul_list.append(image_filename)
            else:
                # Get the image source URL without clicking
                image_src = img.get_attribute('src')
                image_filename = os.path.join(current_patent_directory, f"{patent_name}_image_{index}.png")
                download_image(image_src, image_filename)
                # print(f"Image saved: {image_filename}")
                structural_formul_list.append(image_filename)
        except Exception as e:
            print("Failed to save image.")
            print(e)
            continue
    print(structural_formul_list)

    print("--------------")
    print(uipac_list)
    # Close the browser after the loop ends
    browser.quit()


    # Занесение в базу данных
    # Подключаемся к базе данных
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    # Патент, с которым связываем соединения
    patent_number = patent_name
    # Получаем patent_id по patent_number
    patent_id = get_patent_id(cursor, patent_number)




    # Заносим формулы UIPAC в таблицу compounds
    for uipac_formula in uipac_list:
        cursor.execute("INSERT INTO compounds (uipac) VALUES (?)", (uipac_formula,))
        compound_id = cursor.lastrowid
        # Связываем соединение с патентом
        cursor.execute("INSERT INTO compoundsInPatent (compound_id, patent_id) VALUES (?, ?)", (compound_id, patent_id))

    # Заносим пути к изображениям в таблицу compounds
    # Добавляем новые соединения в таблицу compounds
    for structural_formula_path in structural_formul_list:
        cursor.execute("INSERT INTO compounds (structural_formula) VALUES (?)", (structural_formula_path,))
        compound_id = cursor.lastrowid
        cursor.execute("INSERT INTO compoundsInPatent (compound_id, patent_id) VALUES (?, ?)", (compound_id, patent_id))

    # Коммитим изменения и закрываем соединение
    conn.commit()
    conn.close()





# parse_patent("US20230136730A1")

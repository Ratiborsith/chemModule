import os
import pyautogui
import time
time.sleep(3)  # Wait for the image to load
import requests
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
import urllib.request
import pandas as pd
import os
import subprocess
import sys
import pkg_resources
from selenium.common.exceptions import TimeoutException


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


url = "https://patents.google.com/"
chrome_options = webdriver.ChromeOptions()
chrome_options.add_argument("--headless")  # Runs Chrome in headless mode.
chrome_options.add_argument("--disable-gpu")  # Disables GPU hardware acceleration. If software renderer is not in place, then the headless browser will not launch on Windows.
chrome_options.add_argument("--no-sandbox")  # Bypass OS security model. This can be required if running as root on Linux.
chrome_options.add_argument("--disable-dev-shm-usage")  # Overcome limited resource problems.

service = Service()
browser = webdriver.Chrome(service=service, options=chrome_options)
browser.get(url)
df = pd.read_excel("C:\pycharmpro\chemModule\All_patent_all.xlsx")
patent_numbers = df['patent_number'].tolist()

main_directory = "Patent_images"
if not os.path.exists(main_directory):
    os.makedirs(main_directory)

patents_with_no_images = []

# Assuming 'patent_numbers' is a list of patent numbers and 'url' is the URL to start with
for patent in patent_numbers:
    browser.get(url)
    # Wait for the search box to be present and interact with it
    search_box = WebDriverWait(browser, 10).until(
        EC.presence_of_element_located((By.ID, "searchInput")))
    search_box.clear()
    search_box.send_keys(patent)
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
            print(f"No images found for patent {patent}.")
            patents_with_no_images.append(patent)
            continue

    print(f"Total images found for patent {patent}: {len(images)}")
    current_patent_directory = os.path.join(main_directory, patent)

    if not os.path.exists(current_patent_directory):
        os.makedirs(current_patent_directory)

for index, img in enumerate(images, start=1):
    try:
        if flag == 1:
            image_src = img.get_attribute('src')
            image_filename = os.path.join(current_patent_directory, f"{patent}_image_{index}.png")
            download_image(image_src, image_filename)
            #print(f"Image saved: {image_filename}")
        else:
            # Get the image source URL without clicking
            image_src = img.get_attribute('src')
            image_filename = os.path.join(current_patent_directory, f"{patent}_image_{index}.png")
            download_image(image_src, image_filename)
            #print(f"Image saved: {image_filename}")
    except Exception as e:
        print("Failed to save image.")
        print(e)
        continue

# Close the browser after the loop ends
browser.quit()
if patents_with_no_images:
    print("Patents with no images found:", patents_with_no_images)

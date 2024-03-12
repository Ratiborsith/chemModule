import base64
import os
import sqlite3
from io import BytesIO
from flask import Flask, render_template, request, redirect, url_for, jsonify
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, FileField
from wtforms.validators import DataRequired
from wtforms import SelectField
import cirpy    # библиотека для обработки uipac

import time
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
import urllib.request
from selenium.common.exceptions import TimeoutException
import requests
import re
from bs4 import BeautifulSoup
from itertools import product   # для перебора

# -----------------------------------------------------
# функции для сравнения патентов


def calculate_tanimoto_similarity_Smiles(molecule1, molecule2):
    """
    Calculates the Tanimoto similarity coefficient between two molecules.
    """
    # Convert SMILES strings to RDKit molecules
    mol1 = Chem.MolFromSmiles(molecule1)
    mol2 = Chem.MolFromSmiles(molecule2)

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid input provided.")

    # Generate Morgan fingerprints for the molecules
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)

    # Calculate Tanimoto similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

    return similarity


def print_compound_info(connection, compound_id):
    compounds_in_list = []
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM compounds WHERE id=?", (compound_id,))
    compound = cursor.fetchone()
    if compound:
        # мы должны любое представление преобразовать в smiles и добавить в список уже в таком виде
        if compound[1] != 'not':
            print(f"Compound ID: {compound[0]}")
            print(f"SmilesName: {compound[1]}")
            compounds_in_list.append(compound[1])

        if compound[2] != 'not':
            print(f"Compound ID: {compound[0]}")
            print(f"inChiName: {compound[2]}")

            # преобразуем в Smiles
            compoundSmilesNameMol = Chem.MolFromInchi(compound[2])
            compoundSmilesName = Chem.MolToSmiles(compoundSmilesNameMol)
            compounds_in_list.append(compoundSmilesName)

        if compound[3] != 'not':
            print(f"Compound ID: {compound[0]}")
            print(f"molfile_path: {compound[3]}")

            # преобразуем в Smiles

            MOL1 = compound[3]
            # формируем молекулу из МОЛОВ
            molfiles_dir = os.path.join(os.getcwd())
            MOL1 = MOL1.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL1
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)
            if smiles_from_molfile != None:
                compounds_in_list.append(smiles_from_molfile)

        if compound[4] != 'not':
            print(f"Compound ID: {compound[0]}")
            print(f"uipac: {compound[4]}")

            # преобразуем uipac в smiles
            UIPAC = compound[4]
            try:
                smilesFromUIPAC = cirpy.resolve(UIPAC, "smiles")
                if smilesFromUIPAC is not None:
                    compounds_in_list.append(smilesFromUIPAC)
            except:
                pass

    return compounds_in_list

# вспомогательная функция для функции сравнения патентов
def get_compound_formulas(connection, patent_id):
    all_compounds_in_list = []
    cursor = connection.cursor()
    cursor.execute("SELECT compound_id FROM compoundsInPatent WHERE patent_id=?", (patent_id,))
    compound_ids = cursor.fetchall()
    if compound_ids:
        for compound_id in compound_ids:
            all_compounds_in_list += print_compound_info(connection, compound_id[0])
    else:
        print("No compounds found for patent_id:", patent_id)

    return all_compounds_in_list

def get_similarity_patents(patent1_id, patent2_id, comparison_type):


    connection = sqlite3.connect('patents.db')  # Замените на имя вашей базы данных
    all_compounds_in_patent1 = get_compound_formulas(connection, patent1_id)
    all_compounds_in_patent2 = get_compound_formulas(connection, patent2_id)
    connection.close()
    print(all_compounds_in_patent1)
    print(all_compounds_in_patent2)


    # Создаем список пар элементов из обоих списков
    pairs = list(product(all_compounds_in_patent1, all_compounds_in_patent2))

    # Создаем пустой словарь для хранения сходства между парами
    similarity_scores = {}

    # Проходим по каждой паре и вычисляем сходство
    for pair in pairs:
        similarity = calculate_tanimoto_similarity_Smiles(pair[0], pair[1])
        similarity_scores[(pair[0], pair[1])] = similarity

    # Сортируем словарь по значению сходства в порядке убывания
    sorted_similarity = sorted(similarity_scores.items(), key=lambda x: x[1], reverse=True)

    # строчка для вывода сходства на страницу
    SimilarityString = ""

    if comparison_type == "IUPAC":
        # Выводим результаты
        for pair, similarity in sorted_similarity:
            # Для первого соединения
            first_compound = cirpy.resolve(pair[0], "iupac_name")
            # если не нашлось, то оставляем в smiles
            if first_compound == None:
                first_compound = pair[0]

            # Для второго соединения
            # преобразуем в iupac
            second_compound = cirpy.resolve(pair[1], "iupac_name")
            # если не нашлось, то оставляем в smiles
            if second_compound == None:
                second_compound = pair[1]
            SimilarityString += "<br>"
            SimilarityString += "<br>"
            SimilarityString += f"Сходство между {first_compound} и {second_compound}: {similarity}"
            # Печатаем сходство
            print(f"Сходство между {first_compound} и {second_compound}: {similarity}")

    # возвращаем строчку, содержащую весь нужный нам вывод
    elif comparison_type == "SMILES":
        # Выводим результаты
        for pair, similarity in sorted_similarity:
            first_compound = pair[0]
            second_compound = pair[1]

            SimilarityString += "<br>"
            SimilarityString += "<br>"
            SimilarityString += f"Сходство между {first_compound} и {second_compound}: {similarity}"
            # Печатаем сходство
            print(f"Сходство между {first_compound} и {second_compound}: {similarity}")
    return SimilarityString

# конец функций для сравнения патентов
# -----------------------------------------------------

# функции для парсинга
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

# Главная функция для парсинга патентов
def parse_patent(patent_number_url):

    output = "" # строка вывода

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
                output += type_name + ":"
                for formula in chemical_info[type_name]:
                    try:
                        print(formula)
                        output += formula + "<br>"
                        uipac_list.append(formula)
                        i = i + 1
                    except:
                        pass

        prStr = f"Количество найденных формул: {i}"
        print(prStr)
        #output += prStr + "<br>"


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
            prStr = f"В патенте {patent_name} не найдено структурных представлений."
            print(prStr)
            output += "<br>"+prStr + "<br>"
    try:
        prStr = f"Всего для патента {patent_name} найдено изображений: {len(images)}"
        print(prStr)
        #output += "<br>" + prStr + "<br>"
    except:
        pass

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
                output += "<br>" + f"Image saved: {image_filename}" + "<br>"
                structural_formul_list.append(image_filename)
            else:
                # Get the image source URL without clicking
                image_src = img.get_attribute('src')
                image_filename = os.path.join(current_patent_directory, f"{patent_name}_image_{index}.png")
                download_image(image_src, image_filename)
                # print(f"Image saved: {image_filename}")
                output += "<br>" + f"Image saved: {image_filename}" + "<br>"
                structural_formul_list.append(image_filename)
        except Exception as e:
            print("Failed to save image.")
            output += "<br>" + "Failed to save image." + "<br>"
            print(e)
            continue

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


    return output

#import ParseAllFormuls  # Импортируем функции парсинга

# Функция для получения списка патентов из базы данных
def get_patents():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('SELECT id, name, PatentNumber, link FROM patents')
    patents = cursor.fetchall()
    conn.close()
    return patents

def get_patents_MOL():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('SELECT id, name, PatentNumber, link FROM patents')
    patents = cursor.fetchall()
    conn.close()
    return patents

# Определение функции, которая конвертирует изображение в строку base64
def img_to_base64(img):
    img_buffer = BytesIO()
    img.save(img_buffer, format="PNG")
    img_str = base64.b64encode(img_buffer.getvalue()).decode("utf-8")
    return img_str

# Определение функции, которая конвертирует molfile в SMILES
def molfile_to_smiles(molfile_content):
    mol = Chem.MolFromMolBlock(molfile_content)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        return smiles
    return None

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your_secret_key'
login_manager = LoginManager()
login_manager.init_app(app)

class User(UserMixin):
    pass

@login_manager.unauthorized_handler
def unauthorized_callback():
    return redirect(url_for('login'))

@login_manager.user_loader
def load_user(user_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM users WHERE id = ?', (user_id,))
    user_data = cursor.fetchone()
    conn.close()

    if user_data:
        user = User()
        user.id = user_data[0]  # Предполагается, что первый столбец - это ID пользователя
        return user
    else:
        return None

# форма сравнения
class ChemComparisonForm(FlaskForm):
    patent_id = SelectField('Select Patent', coerce=int, validators=[DataRequired()])
    chem_input2 = StringField('SMILES 2', validators=[DataRequired()])
    submit = SubmitField('Compare')

# форма парсинга
class ChemParseForm(FlaskForm):
    patent_id = SelectField('Select Patent', coerce=int, validators=[DataRequired()])
    submit = SubmitField('Parse')


def calculate_tanimoto_similarity(molecule1, molecule2):
    """
    Calculates the Tanimoto similarity coefficient between two molecules.
    """
    # Convert SMILES strings to RDKit molecules
    mol1 = Chem.MolFromMolBlock(molecule1)
    mol2 = Chem.MolFromSmiles(molecule2)

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid input provided.")

    # Generate Morgan fingerprints for the molecules
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)

    # Calculate Tanimoto similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

    return similarity

# Маршрут для отображения списка патентов
@app.route('/')
def mainPage():
    patents = get_patents()
    return render_template('main.html', patents=patents)


@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']

        conn = sqlite3.connect('patents.db')
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM users WHERE login = ? AND password = ?', (username, password))
        user_data = cursor.fetchone()
        conn.close()

        if user_data:
            user = User()
            user.id = user_data[0]  # Предполагается, что первый столбец - это ID пользователя
            login_user(user)
            return redirect(url_for('index'))  # Перенаправляем на страницу '/similarity'
        else:
            error = "Неверный логин или пароль"

    return render_template('login.html', error=error)


@app.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))


# функция с основными рассчетами
@app.route('/process_comparison', methods=['POST'])
@login_required
def process_comparison():
    form = ChemComparisonForm(request.form)
    similarity = None
    img_base64_1 = None
    img_base64_2 = None
    smiles_from_molfile = None
    selected_second_patent_id = request.form.get('second_patent_id_hidden')

    # Получаем значения из скрытых полей формы

    # номера текущих слайдов
    slick1 = int(request.form.get('slick-active1'))  # на нулевом слайде MOL на единичном слайде SMILES ВСЕГДА на обоих слайдерах
    slick2 = int(request.form.get('slick-active2'))     # на 3 слайде InChi

    # значения изнутри слайдов

    # MOL для 1-го
    MOL1 = request.form.get('molfile_path_hidden')

    # SMILES для 1-го
    SMILES1 = request.form.get('chem_input1_hidden')

    # INCHI для 1-го
    INCHI1 = request.form.get('chem_input_InChi1_hidden')

    # UIPAC для 1-го
    UIPAC1 = request.form.get('chem_input_UIPAC1_hidden')

    # SMILES для 2-го
    SMILES2 = request.form.get('chem_input2_hidden')

    # MOL для 2-го
    MOL2 = request.form.get('molfile_path2_hidden')

    # INCHI для 2-го
    INCHI2 = request.form.get('chem_input_InChi2_hidden')

    # UIPAC для 2-го
    UIPAC2 = request.form.get('chem_input_UIPAC2_hidden')

    if slick1 == 0 and slick2 == 0:
        # сравнение MOL1 и MOL2
        try:
            # формируем первую молекулу из МОЛОВ
            molfiles_dir = os.path.join(app.root_path)
            MOL1 = MOL1.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL1
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)


            # формируем вторую молекулу из МОЛОВ
            molfiles_dir = os.path.join(app.root_path)
            MOL2 = MOL2.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL2
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol2 = Chem.MolFromMolBlock(molfile_content)


            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})


    elif slick1 == 0 and slick2 == 1:
        # сравнение MOL1 и SMILES2
        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL1 = MOL1.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL1
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            mol2 = Chem.MolFromSmiles(SMILES2)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})

    elif slick1 == 1 and slick2 == 0:
        # сравнение SMILES1 и MOL2
        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL2 = MOL2.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL2
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol2 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            mol1 = Chem.MolFromSmiles(SMILES1)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})

    elif slick1 == 1 and slick2 == 1:
        # сравнение SMILES1 и SMILES2
        mol1 = Chem.MolFromSmiles(SMILES1)
        smiles_from_molfile = SMILES1   # ПОТОМ УБРАТЬ

        mol2 = Chem.MolFromSmiles(SMILES2)

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)

    elif slick1 == 2 and slick2 == 1:
        # INCHI1 и SMILES2
        mol1 = Chem.MolFromInchi(INCHI1)
        mol2 = Chem.MolFromSmiles(SMILES2)
        smiles_from_molfile = "NOTHING"

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)


    elif slick1 == 2 and slick2 == 0:
        # INCHI1 и MOL2
        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL2 = MOL2.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL2
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol2 = Chem.MolFromMolBlock(molfile_content)

            mol1 = Chem.MolFromInchi(INCHI1)
            smiles_from_molfile = "NOTHING"

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})



    elif slick1 == 1 and slick2 == 2:
        # SMILES1 и INCHI2
        mol1 = Chem.MolFromSmiles(SMILES1)
        mol2 = Chem.MolFromInchi(INCHI2)
        smiles_from_molfile = "NOTHING"

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)



    elif slick1 == 0 and slick2 == 2:
        # MOL1 и INCHI2
        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL1 = MOL1.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL1
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            mol2 = Chem.MolFromInchi(INCHI2)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})

    elif slick1 == 2 and slick2 == 2:
        # InChi и InChi
        mol1 = Chem.MolFromInchi(INCHI1)
        mol2 = Chem.MolFromInchi(INCHI2)
        smiles_from_molfile = "NOTHING"


        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)


    # всё с uipac
    elif slick1 == 0 and slick2 == 3:
        # MOL1 and UIPAC2
        mol2_SMILES = cirpy.resolve(UIPAC2, 'smiles')

        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL1 = MOL1.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL1
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            mol2 = Chem.MolFromSmiles(mol2_SMILES)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})


    elif slick1 == 1 and slick2 == 3:
        # SMILES1 and UIPAC2
        mol2_SMILES = cirpy.resolve(UIPAC2, 'smiles')

        mol1 = Chem.MolFromSmiles(SMILES1)
        smiles_from_molfile = SMILES1   # ПОТОМ УБРАТЬ

        mol2 = Chem.MolFromSmiles(mol2_SMILES)

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)



    elif slick1 == 3 and slick2 == 0:
        # UIPAC1 and MOL2
        mol1_SMILES = cirpy.resolve(UIPAC1, 'smiles')

        try:
            molfiles_dir = os.path.join(app.root_path)
            MOL2 = MOL2.replace('/', os.sep)
            selected_molfile_abs_path = molfiles_dir + MOL2
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data

            mol2 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            mol1 = Chem.MolFromSmiles(mol1_SMILES)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return jsonify({'error': str(e)})

    elif slick1 == 3 and slick2 == 1:
        # UIPAC1 and SMILES2
        mol1_SMILES = cirpy.resolve(UIPAC1, 'smiles')

        mol1 = Chem.MolFromSmiles(mol1_SMILES)
        smiles_from_molfile = 'NOTHING'   # ПОТОМ УБРАТЬ

        mol2 = Chem.MolFromSmiles(SMILES2)

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)

    elif slick1 == 3 and slick2 == 3:
        # UIPAC1 и UIPAC2
        # преобразуем uipac в smiles
        mol1_SMILES = cirpy.resolve(UIPAC1, 'smiles')
        mol2_SMILES = cirpy.resolve(UIPAC2, 'smiles')

        # дальше проводим обычные операции, как раньше над smiles

        # сравнение SMILES1 и SMILES2
        mol1 = Chem.MolFromSmiles(mol1_SMILES)
        smiles_from_molfile = mol1_SMILES   # ПОТОМ УБРАТЬ

        mol2 = Chem.MolFromSmiles(mol2_SMILES)

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)

    elif slick1 == 2 and slick2 == 3:
        # inchi1 and uipac2
        mol2_SMILES = cirpy.resolve(UIPAC2, 'smiles')
        mol1 = Chem.MolFromInchi(INCHI1)
        mol2 = Chem.MolFromSmiles(mol2_SMILES)
        smiles_from_molfile = "NOTHING"

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)


    elif slick1 == 3 and slick2 == 2:
        # uipac1 and inchi2
        mol1_SMILES = cirpy.resolve(UIPAC1, 'smiles')
        mol1 = Chem.MolFromSmiles(mol1_SMILES)
        mol2 = Chem.MolFromInchi(INCHI2)
        smiles_from_molfile = "NOTHING"

        if mol1 is not None and mol2 is not None:
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

            img1 = Draw.MolToImage(mol1)
            img_base64_1 = img_to_base64(img1)

            img2 = Draw.MolToImage(mol2)
            img_base64_2 = img_to_base64(img2)

    return jsonify({'similarity': similarity, 'img1': img_base64_1, 'img2': img_base64_2, 'smiles_from_molfile': smiles_from_molfile})


# страница управления патентами
@app.route('/parse', methods=['GET', 'POST'])
@login_required
def parse():
    form = ChemParseForm()
    patents = get_patents()
    form.patent_id.choices = [(patent[0], f"{patent[1]} - {patent[2]}") for patent in patents]
    return render_template('parse.html', form=form)

# обработчик для парсинга патента

@app.route('/parse_patent', methods=['POST'])
def parse_patentPage():
    patent_id = request.json.get('patent_id')  # Используйте request.json для получения данных JSON
    user_id = 1

    # Устанавливаем соединение с базой данных
    connection = sqlite3.connect('patents.db')
    cursor = connection.cursor()

    # Выполняем SQL-запрос для получения PatentNumber по id
    cursor.execute("SELECT PatentNumber FROM patents WHERE id=?", (patent_id,))
    PatentNumber = cursor.fetchone()

    # Закрываем соединение
    connection.close()

    parsed_data = parse_patent(PatentNumber[0])
    parsed_data = parsed_data.replace('\n', '<br>')


    return jsonify(parsed_data)  # Верните результаты парсинга как JSON


# ----------------------------
# страница сравнения патентов
@app.route('/patentsSimilarity', methods=['GET', 'POST'])
@login_required
def patentsSimilarity():
    form = ChemParseForm()
    patents = get_patents()
    form.patent_id.choices = [(patent[0], f"{patent[1]} - {patent[2]}") for patent in patents]
    return render_template('patentsSimilarity.html', form=form)


# обработчик для сравнения двух патентов

@app.route('/patents_similar', methods=['POST'])
def patents_similar():
    patent1_id = request.json.get('patent1_id')  # Получим id 1-го патента
    patent2_id = request.json.get('patent2_id')  # Получим id 2-го патента
    comparison_type = request.json.get('comparison_type')     # получим значение радиобатонов


    user_id = 1

    # Устанавливаем соединение с базой данных
    connection = sqlite3.connect('patents.db')
    cursor = connection.cursor()

    # Выполняем SQL-запрос для получения PatentNumber по id
    cursor.execute("SELECT PatentNumber FROM patents WHERE id=?", (patent1_id,))
    Patent1Number = cursor.fetchone()

    # для второго патента
    cursor.execute("SELECT PatentNumber FROM patents WHERE id=?", (patent2_id,))
    Patent2Number = cursor.fetchone()
    # Закрываем соединение
    connection.close()

    comparison_results = f"Результат сравнения двух патентов {Patent1Number[0]} и {Patent2Number[0]}:"

    comparison_results += get_similarity_patents(patent1_id, patent2_id, comparison_type)    # вычисление схожести формул в патенте

    return jsonify(comparison_results)  # Верните результаты парсинга как JSON


# страница сравнения химических соединений
@app.route('/similarity', methods=['GET', 'POST'])
@login_required
def index():
    form = ChemComparisonForm()
    patents = get_patents()
    form.patent_id.choices = [(patent[0], f"{patent[1]} - {patent[2]}") for patent in patents]
    return render_template('index.html', form=form)


@app.route('/get_molfiles_for_patent/<int:patent_id>')
@login_required
def get_molfiles_for_patent(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('SELECT molfile_path FROM compounds WHERE id IN (SELECT compound_id FROM compoundsInPatent WHERE patent_id = ?) AND molfile_path != "not"', (patent_id,))
    molfiles_for_patent = cursor.fetchall()

    conn.close()

    return jsonify({'molfiles': molfiles_for_patent})

@app.route('/get_smiles_for_patent/<int:patent_id>')
@login_required
def get_smiles_for_patent(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('SELECT smilesName FROM compounds WHERE id IN (SELECT compound_id FROM compoundsInPatent WHERE patent_id = ?) AND smilesName != "not"', (patent_id,))
    smiles_for_patent = cursor.fetchall()

    conn.close()

    return jsonify({'smiles': [smile[0] for smile in smiles_for_patent]})

@app.route('/get_INCHI_for_patent/<int:patent_id>')
@login_required
def get_INCHI_for_patent(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('SELECT inChiName FROM compounds WHERE id IN (SELECT compound_id FROM compoundsInPatent WHERE patent_id = ?) AND inChiName != "not"', (patent_id,))
    inchi_for_patent = cursor.fetchall()

    conn.close()

    return jsonify({'inchies': inchi_for_patent})

# UIPAC получение соединений в список
@app.route('/get_UIPAC_for_patent/<int:patent_id>')
@login_required
def get_UIPAC_for_patent(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('SELECT uipac FROM compounds WHERE id IN (SELECT compound_id FROM compoundsInPatent WHERE patent_id = ?) AND uipac != "not"', (patent_id,))
    UIPAC_for_patent = cursor.fetchall()

    conn.close()

    return jsonify({'UIPACS': UIPAC_for_patent})


@app.route('/get_patents_except/<int:patent_id>')
@login_required
def get_patents_except(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('''SELECT DISTINCT p.id, p.name, p.link
FROM patents p
JOIN compoundsInPatent cp ON p.id = cp.patent_id
JOIN compounds c ON cp.compound_id = c.id
WHERE c.smilesName != "not";
''')
    patents_except = cursor.fetchall()

    conn.close()

    return jsonify({'patents': patents_except})


@app.route('/get_patents_with_molfiles')
def get_patents_with_molfiles():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT p.id, p.name, p.link
FROM patents p
JOIN compoundsInPatent cp ON p.id = cp.patent_id
JOIN compounds c ON cp.compound_id = c.id
WHERE c.molfile_path != "not"
    ''')
    patents_with_molfiles = cursor.fetchall()
    conn.close()
    return jsonify({'patents': patents_with_molfiles})

@app.route('/get_patents_with_smiles')
def get_patents_with_smiles():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT p.id, p.name, p.link
FROM patents p
JOIN compoundsInPatent cp ON p.id = cp.patent_id
JOIN compounds c ON cp.compound_id = c.id
WHERE c.smilesName != "not"
    ''')
    patents_with_smiles = cursor.fetchall()
    conn.close()
    return jsonify({'patents': patents_with_smiles})

# получить патенты с InChi
@app.route('/get_patents_with_INCHI')
def get_patents_with_INCHI():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT p.id, p.name, p.link
FROM patents p
JOIN compoundsInPatent cp ON p.id = cp.patent_id
JOIN compounds c ON cp.compound_id = c.id
WHERE c.inChiName != "not"
    ''')
    patents_with_inchi= cursor.fetchall()
    conn.close()
    return jsonify({'patents': patents_with_inchi})

# получить патенты с UIPAC
@app.route('/get_patents_with_UIPAC')
def get_patents_with_UIPAC():
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()
    cursor.execute('''
        SELECT DISTINCT p.id, p.name, p.link
FROM patents p
JOIN compoundsInPatent cp ON p.id = cp.patent_id
JOIN compounds c ON cp.compound_id = c.id
WHERE c.uipac != "not"
    ''')
    patents_with_UIPAC= cursor.fetchall()
    conn.close()
    return jsonify({'patents': patents_with_UIPAC})


@app.route('/about')
def about():
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=True)
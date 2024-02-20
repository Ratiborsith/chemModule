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

class ChemComparisonForm(FlaskForm):
    patent_id = SelectField('Select Patent', coerce=int, validators=[DataRequired()])
    chem_input2 = StringField('SMILES 2', validators=[DataRequired()])
    submit = SubmitField('Compare')

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




# страница сравнения химических соединен
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
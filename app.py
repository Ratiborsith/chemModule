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

# Функция для получения списка патентов из базы данных
def get_patents():
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

# страница сравнения химических соединен
@app.route('/similarity', methods=['GET', 'POST'])
@login_required
def index():
    form = ChemComparisonForm()
    similarity = None
    img_base64_1 = None
    img_base64_2 = None
    smiles_from_molfile = None
    molfiles_for_patent = None

    patents = get_patents()
    form.patent_id.choices = [(patent[0], f"{patent[1]} - {patent[2]}") for patent in patents]

    if form.validate_on_submit():
        selected_patent_id = form.patent_id.data
        selected_molfile_path = request.form.get('molfile_path')  # Получаем выбранный пользователем molfile_path
        selected_second_patent_id = request.form.get('second_patent_id')  # Получаем выбр. пользователем второй патент
        chem_input2 = request.form.get('chem_input2')  # Получаем выбранный пользователем SMILES для второго соединения

        try:
            # Получаем путь к каталогу molfiles внутри проекта
            molfiles_dir = os.path.join(app.root_path)
            selected_molfile_path = selected_molfile_path.replace('/', os.sep)
            # Формируем абсолютный путь к выбранному файлу
            selected_molfile_abs_path = molfiles_dir + selected_molfile_path
            with open(selected_molfile_abs_path, 'r') as molfile:
                molfile_data = molfile.read()
                molfile_content = molfile_data  # Считываем содержимое molfile_data

            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)

            # Строим вторую молекулу по введенному пользователем SMILES
            mol2 = Chem.MolFromSmiles(chem_input2)

            if mol1 is not None and mol2 is not None:
                similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), Chem.MolToSmiles(mol2))

                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return render_template('index.html', form=form, error=str(e))

    return render_template('index.html', form=form, similarity=similarity, img1=img_base64_1, img2=img_base64_2,
                           smiles_from_molfile=smiles_from_molfile, molfiles_for_patent=molfiles_for_patent)


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

@app.route('/get_patents_except/<int:patent_id>')
@login_required
def get_patents_except(patent_id):
    conn = sqlite3.connect('patents.db')
    cursor = conn.cursor()

    cursor.execute('SELECT id, name FROM patents WHERE id != ?', (patent_id,))
    patents_except = cursor.fetchall()

    conn.close()

    return jsonify({'patents': patents_except})

@app.route('/about')
def about():
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=True)

import base64
from io import BytesIO
from flask import Flask, render_template, request, redirect, url_for
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, FileField
from wtforms.validators import DataRequired

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
    # Здесь нужно реализовать загрузку пользователя из базы данных, которую я создам в дальнейшем
    # В данном примере загружаем фиктивного пользователя с идентификатором '1'
    user = User()
    user.id = '1'
    return user

class ChemComparisonForm(FlaskForm):
    molfile_input = FileField('Molfile')
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

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        # Неразвитая логика проверки имени пользователя и пароля
        username = request.form['username']
        password = request.form['password']

        # В данном примере допустим, что проверка пройдена
        user = User()
        user.id = '1'
        login_user(user)
        return redirect(url_for('index'))  # Перенаправляем на страницу '/similarity'

    return render_template('login.html')

@app.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('login'))

# страница сравнения химических соединений
@app.route('/similarity', methods=['GET', 'POST'])
@login_required
def index():
    form = ChemComparisonForm()
    similarity = None

    img_base64_1 = None
    img_base64_2 = None
    smiles_from_molfile = None

    if form.validate_on_submit():
        molfile = form.molfile_input.data
        chem_input2 = form.chem_input2.data.strip()

        try:
            molfile_content = molfile.read().decode('utf-8')  # Считываем содержимое molfile
            mol1 = Chem.MolFromMolBlock(molfile_content)
            smiles_from_molfile = molfile_to_smiles(molfile_content)
            similarity = calculate_tanimoto_similarity(Chem.MolToMolBlock(mol1), chem_input2)

            if mol1 is not None:
                img1 = Draw.MolToImage(mol1)
                img_base64_1 = img_to_base64(img1)

            # Convert SMILES strings to RDKit molecules for drawing images
            mol2 = Chem.MolFromSmiles(chem_input2)

            if mol2 is not None:
                img2 = Draw.MolToImage(mol2)
                img_base64_2 = img_to_base64(img2)

        except Exception as e:
            return render_template('index.html', form=form, error=str(e))

    return render_template('index.html', form=form, similarity=similarity, img1=img_base64_1, img2=img_base64_2,
                           smiles_from_molfile=smiles_from_molfile)


@app.route('/about')
def about():
    return render_template('about.html')



if __name__ == '__main__':
    app.run(debug=True)
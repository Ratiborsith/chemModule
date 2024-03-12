import os
import sqlite3
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
import cirpy    # библиотека для обработки uipac
from itertools import product

# Определение функции, которая конвертирует molfile в SMILES
def molfile_to_smiles(molfile_content):
    mol = Chem.MolFromMolBlock(molfile_content)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        return smiles
    return None


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

def get_similarity_patents(patent1_id, patent2_id):


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
        SimilarityString += f"Сходство между {first_compound} и {second_compound}: {similarity}"
        # Печатаем сходство
        print(f"Сходство между {first_compound} и {second_compound}: {similarity}")
    # возвращаем строчку, содержащую весь нужный нам вывод
    return SimilarityString
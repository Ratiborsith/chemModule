from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

import cirpy    # библиотека для обработки uipac
def calculate_tanimoto_similarity(molecule1, molecule2):
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


# представление structural formula
smilesMol1 = str(cirpy.resolve('Methyl 2-(4-(2-(((butylthio)methyl)thio)-3-cyano-6-(thiazol-2-yl)pyridin-4-yl)phenoxy)acetate', 'smiles'))
print(smilesMol1)
print("C(C(=O)O)N")

print(calculate_tanimoto_similarity(smilesMol1, "C(C(=O)O)N"))

# структурная - это графически. Надо просто натренировать модель!!!!!!!!!
# молекулярная не нужна

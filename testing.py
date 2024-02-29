from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs

from glyles import convert

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
    img1 = Draw.MolToImage(mol1)

    img2 = Draw.MolToImage(mol2)

    # Отображение изображений
    img1.show()
    img2.show()

    return similarity


# представление structural formula
IUPAC1MOL1 = '3-[2-(butylsulfanylmethylsulfanyl)-3-cyano-6-(1,3-thiazol-2-yl)pyridin-4-yl]-N,N-dimethylbenzamide'
smilesMol1 = str(cirpy.resolve(IUPAC1MOL1, 'smiles'))
IUPAC2Mol1 = str(cirpy.resolve(smilesMol1, 'iupac_name'))
print(IUPAC1MOL1)
print(IUPAC2Mol1)
print(smilesMol1)
print("C(C(=O)O)N")

#print(calculate_tanimoto_similarity(smilesMol1, "C(C(=O)O)N"))

print(calculate_tanimoto_similarity("CCCS(=O)C1=C(C2=C(C=C(C3=NC=CS3)N=C2S1)C4=CC=C(C=C4)CO)N", "CCCCS(=O)C1=C(C2=C(C=C(C3=NC=CS3)N=C2S1)C4=CC=C(C=C4)CO)N"))



# структурная - это графически. Надо просто натренировать модель!!!!!!!!!
# молекулярная не нужна

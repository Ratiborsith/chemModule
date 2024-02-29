from DECIMER import predict_SMILES
import Pyosra
# Chemical depiction to SMILES translation
image_path = "structuralFormuls/US20230136730A1-20230504-C00205.png"
SMILES = predict_SMILES(image_path)
print(SMILES)
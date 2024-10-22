import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

def sdf_to_csv(sdf_file, csv_file):
    # Cargar el archivo SDF
    sdf_data = Chem.SDMolSupplier(sdf_file)
    
    # Convertir los datos a un DataFrame de pandas
    df = PandasTools.LoadSDF(sdf_file)
    
    # Añadir columna de SMILES
    df['SMILES'] = [Chem.MolToSmiles(mol) if mol else None for mol in sdf_data]
    
    # Guardar a CSV
    df.to_csv(csv_file, sep=';', index=False)
    print(f'Archivo guardado como {csv_file}')

# Ruta del archivo SDF y nombre del archivo CSV de salida
input_sdf = 'C:/Users/Enrique/Documents/New ProtoADME Models/CYP_1A2_substrate/raw_data/Tian database/CypReact_Testing_Set.sdf'
output_csv = 'CypReact_test_dataset.csv'

sdf_to_csv(input_sdf, output_csv)

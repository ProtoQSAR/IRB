############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i)
##############################################################################
#################################### IMPORTS #################################
import os
import re   #"""importar regular expresions"""
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
##############################################################################
################################ INITIAL VARIABLES ###########################

pd.set_option('mode.chained_assignment', None) # pd.set_option('mode.chained_assignment', 'warn')

parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

def sanitize_smiles(smiles):
    """
    Convierte una cadena SMILES a un formato canónico usando RDKit.
    Si el SMILES no es válido, devuelve None.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol)
    else:
        return None
def eliminar_version(model):
    """
    Elimina la versión del modelo en una cadena.
    Ejemplo: 'ModelName (version x.x.x)' -> 'ModelName'
    """
    return re.sub(r"\(version \d*\.\d*\.\d*\)?", "", model).strip()

def detectar_skiprow(file_path, header_start="No.", return_models=False):
    """
    Detecta la línea donde se encuentra el encabezado en un archivo de VEGA
    y obtiene los nombres de los modelos si se requiere.

    Parámetros:
    - file_path: Ruta al archivo de predicciones de VEGA.
    - header_start: Texto que indica el inicio del encabezado (por defecto, 'No.').
    - return_models: Si es True, devuelve la lista de modelos encontrados.

    Retorna:
    - Si return_models es False: índice (int) de la fila donde inicia el encabezado.
    - Si return_models es True: lista de nombres de modelos (list).
    """
    models = []
    with open(file_path, "r", encoding="ANSI") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if header_start in line:
            if return_models:
                # Los modelos están entre la primera línea y dos líneas antes del encabezado
                models = [lines[j].strip() for j in range(1, i - 2)]
                models = [eliminar_version(model) for model in models]
                return models
            return i
    return None if not return_models else models

def read_vega_predictions(file_paths):
    """
    Lee y procesa los archivos de predicciones de VEGA.

    Parámetros:
    - file_paths: Lista de rutas a los archivos de predicciones de VEGA (.txt)

    Retorna:
    - DataFrame de pandas con las predicciones procesadas.
    """
    list_to_merge = []
    for file_path in file_paths:
        # Detectar la fila donde inicia el encabezado y obtener los modelos
        skiprows = detectar_skiprow(file_path)
        models = detectar_skiprow(file_path, return_models=True)

        # Comprobar si se detectó el encabezado
        if skiprows is None:
            print(f"No se encontró el encabezado en el archivo {file_path}.")
            continue

        # Leer el archivo de predicciones
        df = pd.read_csv(
            file_path,
            sep="\t",
            encoding="ANSI",
            skiprows=skiprows,
        )

        # Renombrar columnas de predicciones para identificarlas fácilmente
        rename_columns = {
            df.columns[i]: f'Pred_{df.columns[i]}'
            for i in range(4, len(df.columns), 2)
        }
        df.rename(columns=rename_columns, inplace=True)

        new_columns_exp = {}
        new_columns_AD = {}

        # Procesar columnas experimentales y de AD
        for i in range(3, len(df.columns), 2):
            col_name = df.columns[i]
            # Procesar valores experimentales
            new_columns_exp[f'Exp_{col_name}'] = np.where(
                df[col_name].str.contains("EXPERIMENTAL", na=False),
                df[col_name].str.extract(r'(-?[\d.]+)')[0].astype(float),
                np.nan
            )
            # Procesar AD: 0 si contiene "(LOW reliability)", 1 en caso contrario
            new_columns_AD[f'AD_{col_name}'] = np.where(
                df[col_name].str.contains(r'\(LOW reliability\)', na=False),
                0,
                1
            )

        # Añadir las nuevas columnas al DataFrame
        for col_name, col_data in new_columns_exp.items():
            df[col_name] = col_data

        for col_name, col_data in new_columns_AD.items():
            df[col_name] = col_data

        # Filtrar las columnas relevantes
        filtered_columns = df.filter(regex=r'^(tId|SMILES|Exp_|AD_|Pred_)', axis=1)
        df = filtered_columns

        # Añadir el DataFrame procesado a la lista
        list_to_merge.append(df)

    # Concatenar todos los DataFrames
    if list_to_merge:
        predictions_df = pd.concat(list_to_merge, ignore_index=True)
    else:
        predictions_df = pd.DataFrame()
        print("No se pudieron procesar archivos de VEGA.")
    
    return predictions_df
def read_protopred_predictions(file_paths):
    """
    Lee y procesa los archivos de predicciones de ProtoPRED.

    Parámetros:
    - file_paths: Lista de rutas a los archivos CSV de predicciones de ProtoPRED.

    Retorna:
    - DataFrame de pandas con las predicciones procesadas.
    """
    list_to_merge = []
    for file_path in file_paths:
        for file_path in file_paths:
        # Leer el archivo Excel
        # Suponiendo que el nombre de la hoja es 'Results' o puede ser el nombre del endpoint específico
            try:
                df = pd.read_excel(file_path, sheet_name='Water solubility')
            except ValueError:
                # Si la hoja 'Results' no existe, leer la primera hoja
                df = pd.read_excel(file_path)
        
        # Renombrar columnas para estandarizar
        df.rename(columns={
            'Other Regulatory ID': 'ID_ProtoPRED',
            'SMILES': 'SMILES',
            'Predicted numerical (model units)': 'ProtoPRED_pred',
            'Applicability domain**': 'Applicability_domain',
            # Añade más columnas si es necesario
        }, inplace=True)
        print(df)
        # Procesar el Applicability Domain (AD)
        if 'Applicability_domain' in df.columns:
            df['ProtoPRED_AD'] = df['Applicability_domain'].apply(
                lambda x: 1 if 'Inside' in str(x) else 0
            )
            df.drop(columns=['Applicability_domain'], inplace=True)

        # Procesar valores de predicción
        # Si necesitas hacer alguna conversión adicional, puedes agregarla aquí

        # Sanear los SMILES
        df['SMILES'] = df['SMILES'].apply(sanitize_smiles)

        # Eliminar duplicados basados en 'SMILES'
        df.drop_duplicates(subset='SMILES', inplace=True)

        # Añadir el DataFrame procesado a la lista
        list_to_merge.append(df)

    # Concatenar todos los DataFrames
    if list_to_merge:
        predictions_df = pd.concat(list_to_merge, ignore_index=True)
    else:
        predictions_df = pd.DataFrame()
        print("No se pudieron procesar archivos de ProtoPRED.")

    return predictions_df

if __name__ == "__main__":
    import glob

    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Procesar archivos de VEGA
    vega_files = glob.glob(os.path.join(current_dir, '*.txt'))
    if not vega_files:
        print("No se encontraron archivos de predicciones de VEGA en el directorio actual.")
    else:
        vega_predictions_df = read_vega_predictions(vega_files)
        # Guardar o manejar el DataFrame de VEGA según lo necesites

    # Procesar archivos de ProtoPRED
    protopred_files = glob.glob(os.path.join(current_dir, 'protopred_*.xlsx'))
    if not protopred_files:
        print("No se encontraron archivos de predicciones de ProtoPRED en el directorio actual.")
    else:
        protopred_predictions_df = read_protopred_predictions(protopred_files)
        # Guardar o manejar el DataFrame de ProtoPRED según lo necesites

    # Fusionar los DataFrames si ambos existen
    if not vega_predictions_df.empty and not protopred_predictions_df.empty:
        merged_df = pd.merge(vega_predictions_df, protopred_predictions_df, on='SMILES', how='inner')
        output_file = os.path.join(current_dir, 'Merged_predictions.csv')
        print(merged_df.shape)
        # merged_df.to_csv(output_file, sep=';', index=False)
        print(f"Predicciones fusionadas guardadas en {output_file}")
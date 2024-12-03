#### First we import the required modules and functions

import pandas as pd
import json
import numpy as np
from rdkit import Chem
import re

## Create CONSTANTS
QSAR_MODELS = {
    "ProtoPRED" : ["melting_point",
                "boiling_point",
                "vapour_pressure",
                "surface_tension",
                "water_solubility",
                "log_kow",
                "log_d",
                "skin_irritation",
                "eye_irritation",
                "skin_sensitisation",
                "mutagenicity_ames",
                "chromosomal_aberration",
                "hprt_assay",
                "genotoxicity_micronucleus",
                "comet_assay",
                "acute_oral",
                "developmental_toxicity",
                "carcinogenicity",
                "neurotoxicity",
                "ready_biodegradability",
                "fish_acute_toxicity",
                "sludge_inhibition",
                "persistence_water",
                "persistence_soil",
                "persistence_sediment",
                "sorption",
                "bioconcentration_factor",

                "bioavailability20",
                "bioavailability30",
                "caco-2_permeability",
                "p-gp_inhibitor",
                "p-gp_substrate",
                "skin_permeability",
                "human_intestinal_absorption",

                "CYP450_1A2_inhibitor",
                "CYP450_1A2_substrate",
                "CYP450_2C19_inhibitor",
                "CYP450_2C19_substrate",
                "CYP450_2C9_inhibitor",
                "CYP450_2C9_substrate",
                "CYP450_2D6_inhibitor",
                "CYP450_2D6_substrate",
                "CYP450_3A4_inhibitor",
                "CYP450_3A4_substrate",

                "blood-brain_barrier",
                "plasma-protein_binding",
                "volume_of_distribution",

                "half-life",
                "human_liver_microsomal",
                "OATP1B1",
                "OATP1B3",
                "BSEP"],
    "ProtoAquaTox" : ["Acute toxicity on invertebrates",
                      "Acute growth inhibition on algae",
                      "Long-term toxicity on fish",
                      "Long-term toxicity on invertebrates",
                      "Chronic growth inhibition on algae"],
    "Vega" : ["Earthworm toxicity",
              "Plasma Protein Binding",
              "Bee acute toxicity"],
    "TEST" : ["Flash_point",
              "Viscosity_at_25C"],
    "KATE" : ["fels"],
    "MLTOX" : ["phototoxicity"],
    "EPISUITE" :["OH_rate_constant",
                 "Ozone_rate_estimation"],
    "QSARToolbox" : ["ER activation", "ER binding"],
    "ADMETlab" :["ROA", "Respiratory"],
    "STopTox" : ["acute_dermal"]
}

def main():
    global OUTPUTS_PATH
    global ProtoTOX_file
    global ProtoECO_file
    global ProtoPHYSCHEM_file
    global ProtoADME_file
    global ProtoAquaTox_file
    global VEGA_file
    global MLTOX_file
    global MLTOX_file
    global EPISUITE_file
    global QTB_file
    global ADMETlab_file
    OUTPUTS_PATH = input("Introduce the PATH of the output files (enter to: ../data/outputs/): ") or "../data/outputs/"
    ProtoTOX_file = input("Introduce the name of the ProtoTOX file (enter to: ProtoTOX_output.xlsx): ") or "ProtoTOX_output.xlsx"
    ProtoECO_file = input("Introduce the name of the ProtoECO file (enter to: ProtoECO_output.xlsx): ") or "ProtoECO_output.xlsx"
    ProtoPHYSCHEM_file = input("Introduce the name of the ProtoPHYSCHEM file (enter to: ProtoPHYSCHEM_output.xlsx): ") or "ProtoPHYSCHEM_output.xlsx"
    ProtoADME_file = input("Introduce the name of the ProtoADME file (enter to: ProtoADME_output.xlsx): ") or "ProtoADME_output.xlsx"
    ProtoAquaTox_file = input("Introduce the name of the ProtoAquaTox file (enter to: input_ProtoAquaTox_UserDefinedTestDetailedResults.xlsx): ") or "input_ProtoAquaTox_UserDefinedTestDetailedResults.xlsx"
    VEGA_file = input("Introduce the name of the VEGA file (enter to: report_summary.txt): ") or "report_summary.txt"
    for endpoint in QSAR_MODELS['TEST']:
        exec(f"{endpoint}_file = input('Introduce the name of the {endpoint} file (enter to: {endpoint}_Consensus.csv): ') or endpoint+'_Consensus.csv'")
    MLTOX_file = input("Introduce the name of the MLTOX file (enter to: MLTox_output.csv): ") or "MLTox_output.csv"
    EPISUITE_file = input("Introduce the name of the EPISUITE file (enter to: BATCH01.OUT): ") or "BATCH01.OUT"
    QTB_file = input("Introduce the name of the QSARToolbox file (enter to: QTB_estrogen.csv): ") or "QTB_estrogen.csv"
    ADMETlab_file = input("Introduce the name of the ADMETlab file (enter to: ADMETLab_output.csv): ") or "ADMETLab_output.csv"

    search_predictions_batch()

def igualdad_lista(x):
    lista_de_valores = list(x)
    if all(x == lista_de_valores[0] for x in lista_de_valores):
        return lista_de_valores[0]
    else:
        return 'Inconclusivo'
    
def check_star(value):
    return "Experimental" if "*" in str(value) else "Predicted"

def remove_units_and_stars(value):
    # Verificar si el valor es una cadena
    if isinstance(value, str):
        # Verificar si hay algún dígito en la cadena
        if re.search(r'\d', value):
            # Eliminar unidades y asteriscos, manteniendo solo la parte numérica y el punto decimal
            return ''.join(filter(lambda x: x.isdigit() or x == '.', value))
        else:
            # Si no hay dígitos, se asume que es una variable binaria y se eliminan solo los asteriscos
            return value.replace('*', '')
    return value

def extraer_unidades(df, fila=0):
    '''
    Extrae las unidades de medida de una predicción

    Args:
    - df (DataFrame): DataFrame principal

    - fila (int): fila que se quiere emplear para realizar la comparación

    Returns:
    - unidades  (str): cadena con las unidades
    '''
    # Extra el valor de Predicted value como cadena
    predicted_value_str = str(df.loc[fila, 'Predicted value'])

    # Extra el valor de Predicted without units como cadena
    predicted_without_units_str = str(df.loc[fila, 'Predicted without units'])

    # Compara la longitud de las cadena, la de mayor longitud contiene unidades
    if len(predicted_value_str) > len(predicted_without_units_str):
        # Elimina valores, ya que los tienen en común y los espacios en blanco
        unidades = predicted_value_str.replace(predicted_without_units_str,
                                               '').strip()  # Extraer unidades
    else:
        None

    return unidades

def predictions_ProtoPRED(output_ProtoPRED):
    model_endpoints = {
        'Melting point': 'ºC',
        'Boiling point': 'ºC',
        'Vapour pressure': 'mmHg',
        'Water solubility': 'mol/L',
        'Log kow': 'Dimensionless',
        'Log d': 'Dimensionless',
        'Surface tension': 'dyn/cm',
        'Neurotoxicity': 'mg/kg [LD50]',
        'Fish acute toxicity': 'mg/L [LC50]',
        'Persistence water': 'days',
        'Persistence soil': 'days',
        'Persistence sediment': 'days',
        'Sorption': 'L/Kg',
        'Bioconcentration factor': 'L/Kg',
    }
    # Apply the function to all columns except SMILES
    for col in output_ProtoPRED.columns[1:]:
        output_ProtoPRED[f"{col} Type"] = output_ProtoPRED[col].apply(check_star)

    # Apply the function to all columns except SMILES
    new_column_order = []
    for col in output_ProtoPRED.columns:
        if "Type" not in col:  # Skip "Type" columns to avoid duplication
            new_column_order.append(col)
            # Append the "Type" column immediately after the original column, if it exists
            type_col = f"{col} Type"
            if type_col in output_ProtoPRED.columns:
                new_column_order.append(type_col)

    ProtoPRED_df = output_ProtoPRED[new_column_order]

    for col in ProtoPRED_df.columns:
        if "SMILES" not in col and "Type" not in col:
            ProtoPRED_df[col] = ProtoPRED_df[col].apply(remove_units_and_stars)

    ProtoPRED_df.rename(columns={col: f"{col} ({model_endpoints[col]})" for col in model_endpoints if col in ProtoPRED_df.columns}, inplace=True)

    ProtoPRED_df.set_index('SMILES', inplace=True)

    return ProtoPRED_df

def ProtoADME_dict(df):
    '''
    Automatiza la creación de un diccionario a partir de un DataFrame

    Args:
    - df (DataFrame): DataFrame principal que contiene todas las
                           categorías y datos.

    Returns:
    - dict: Diccionario de categorías con las columnas seleccionadas.
    '''
    # Diccionario con las columnas para cada propiedad
    ProtoADME_prop = {}

    # Crea un diccionario con las columnas para cada propiedad
    for keys, dataframe in df.items():
        # Evita las dos primeras hojas del excel
        if keys != 'General information' and keys != 'Druglikeness':
            # Verificar si existe la columna Predicted binary (0/1)
            if 'Predicted binary (0/1)' in dataframe.columns:
                ProtoADME_prop[keys] = ['SMILES',
                                        'Experimental value*',
                                        'Predicted value',
                                        'Predicted binary (0/1)',
                                        'Applicability domain**']

            # Verificar si existe la columna Predicted without units
            elif 'Predicted without units' in dataframe.columns:
                ProtoADME_prop[keys] = ['SMILES',
                                        'Experimental value*',
                                        'Predicted value',
                                        'Predicted without units',
                                        'Applicability domain**']

    return ProtoADME_prop

def predictions_ProtoADME(ADME_dataframe, categories_columns):
    '''
    Automatiza la creación de DataFrames a partir de un DataFrame principal y
    un diccionario. Los nombres de las columnas se actualizan para reflejar la
    categoría a la que pertenecen.

    Args:
    - df (DataFrame): DataFrame principal que contiene todas las
                           categorías y datos.
    - categories_columns (dict): Diccionario que mapea nombres de categorías.

    Returns:
    - dict  (dict): Diccionario de DataFrames, cada uno correspondiente a una
            categoría con las columnas seleccionadas.
    '''

    ADME_df = pd.DataFrame()
    dataframes_dict = {}
    # Crea el nuevo dataframe
    for category, columns in categories_columns.items():
        # Crear el DataFrame y seleccionar columnas
        df = ADME_dataframe[category].loc[:, columns]

        # Inicia la lista de nuevos nombres de columna
        new_column_names = ['Sanitized SMILES']

        # Cambia los nombres a las columnas
        for col in columns[1:]:
            if col == 'Predicted without units':
                # Obtiene unidades
                units = extraer_unidades(df, 0)
                new_name = f'{category} Predicted value ({units})'
                new_column_names.append(new_name)
            elif col == 'Experimental value*':
                new_name = f'{category} Experimental value'
                new_column_names.append(new_name)
            elif col == 'Applicability domain**':
                new_name = f'{category} Applicability domain'
                new_column_names.append(new_name)
            else:
                new_column_name = f"{category} {col}"
                new_column_names.append(new_column_name)

        # Asigna los nuevos nombres al df
        df.columns = new_column_names

        # Eliminate duplicated columns
        eliminate = f"{category} Predicted value"
        df.drop(eliminate, axis=1, inplace=True)

        # Reemplazar espacios por guiones bajos y agregar '_df' al final del
        # nombre de cada DataFrame
        dataframe_name = category.replace(' ', '_').replace('-', '_') + '_df'
        dataframes_dict[dataframe_name] = df

    for key, df in dataframes_dict.items():  # Itera sobre los DataFrames
        # Merge basado en 'Sanitized SMILES'
        df.set_index('Sanitized SMILES', inplace=True)
        ADME_df = ADME_df.join(df, how='outer')

    # ADME_df.rename(columns={'Sanitized SMILES':'SMILES'}, inplace=True)

    # Devuelve el diccionario de dataframes
    return ADME_df

def predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint):
    ENDPOINT_DICT = {
        "Acute toxicity on invertebrates" : {
            'Ef':'Immobilisation', 
            'Tt':'Acute'
        },
        "Acute growth inhibition on algae" : {
            'Ef':'Growth', 
            'Tt':'Acute'
        },
        "Long-term toxicity on fish" : {
            'Ef':'Mortality', 
            'Tt':'Chronic'
        },
        "Long-term toxicity on invertebrates" : {
            'Ef':'Immobilisation', 
            'Tt':'Chronic'
        },
        "Chronic growth inhibition on algae" : {
            'Ef':'Growth', 
            'Tt':'Chronic'
        }
    }
    
    endpoint_df = pd.DataFrame()
    mask_Ef = output_ProtoAquaTox['Ef']==ENDPOINT_DICT[endpoint]['Ef']
    mask_Tt = output_ProtoAquaTox['Tt']==ENDPOINT_DICT[endpoint]['Tt']

    smiles_list = []
    predictions_list = []

    for row in output_ProtoAquaTox[mask_Ef & mask_Tt][['SMILES','FINAL OUTCOME']].iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['SMILES'])))
        predictions_list.append(row[1]['FINAL OUTCOME'])

    endpoint_df['SMILES'] = smiles_list
    endpoint_df[endpoint] = predictions_list

    endpoint_df.set_index('SMILES', inplace=True)

    return endpoint_df

def predictions_Vega(output_Vega):
    VEGA_df = pd.DataFrame()
    smiles_list = []

    for row in output_Vega.iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['SMILES'])))
    VEGA_df['SMILES'] = smiles_list

    for i in range(3,len(output_Vega.columns),2):
        assessment_list = []
        predictions_list = []
        db_pred_column = output_Vega.columns[i].split('-assessment')[0]+' prediction'
        db_assesment_column = output_Vega.columns[i].split('-assessment')[0]+' assessment'
        VEGA_column = output_Vega.columns[i]
        for row in output_Vega.iterrows():
            assessment_list.append(row[1][VEGA_column].split('(')[-1].strip(')'))
            predictions_list.append('('.join(row[1][VEGA_column].split('(')[0:-1]))

        VEGA_df[db_pred_column] = predictions_list
        VEGA_df[db_assesment_column] = assessment_list
    
    VEGA_df.set_index('SMILES', inplace=True)

    return VEGA_df

def predictions_TEST(output_TEST, endpoint):
    TEST_df = pd.DataFrame()
    exp_column = endpoint+output_TEST.columns[5]
    pred_column = endpoint+output_TEST.columns[6]
    smiles_list = []
    predicted_list = []
    experimental_list = []
    for row in output_TEST.iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['Query'])))
        experimental_list.append(row[1][5])
        predicted_list.append(row[1][6])
    TEST_df['SMILES'] = smiles_list
    TEST_df[exp_column] = experimental_list
    TEST_df[pred_column] = predicted_list

    TEST_df.set_index('SMILES', inplace=True)

    return TEST_df

def predictions_MLTOX(output_MLTOX):
    #Esta sección hay que remirarla, ya que Slava hablo de una API para solucionar
    #el problema de que no se pueda meter archivo de input ni obtener archivo de output.
    MLTOX_df = pd.DataFrame()
    smiles_list = []
    predicted_list = []
    for row in output_MLTOX.iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['SMILES'])))
        if row[1]['Phototoxicity']==1:
            predicted_list.append('Phototoxic')
        elif row[1]['Phototoxicity']==0:
            predicted_list.append('Non-phototoxic')
        else:
            predicted_list.append('Non predicted')
    
    MLTOX_df['SMILES'] = smiles_list
    MLTOX_df['Phototoxicity'] = predicted_list

    MLTOX_df.set_index('SMILES', inplace=True)

    return MLTOX_df

def predictions_EPISUITE(output_EPISUITE):
    EPISUITE_df = pd.DataFrame()
    
    output_EPISUITE['SMILES'] = output_EPISUITE['SMILES+Id'].apply(lambda x: x.split()[0])
    output_EPISUITE['Predicted OH Rate Constants (cm3/molecule-sec)'] = output_EPISUITE['OH predicted'].apply(lambda x: x.split('(')[0])
    output_EPISUITE['Experimental OH Rate Constants'] = output_EPISUITE['OH experimental'].apply(lambda x: x.split('(')[0])
    output_EPISUITE['Predicted Ozone Reaction Rate (cm3/molecule-sec)'] = output_EPISUITE['Ozone predicted'].apply(lambda x: x.split('(')[0])
    output_EPISUITE['Experimental Ozone Reaction Rate'] = output_EPISUITE['Ozone experimental'].apply(lambda x: x.split('(')[0])

    smiles_list = []
    OH_predicted = []
    OH_exp = []
    ozone_predicted = []
    ozone_exp = []

    for row in output_EPISUITE.iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['SMILES'])))
        OH_predicted.append(row[1]['Predicted OH Rate Constants (cm3/molecule-sec)'])
        OH_exp.append(row[1]['Experimental OH Rate Constants'])
        ozone_predicted.append(row[1]['Predicted Ozone Reaction Rate (cm3/molecule-sec)'])
        ozone_exp.append(row[1]['Experimental Ozone Reaction Rate'])
    
    EPISUITE_df['SMILES'] = smiles_list
    EPISUITE_df['Predicted OH Rate Constants (cm3/molecule-sec)'] = OH_predicted
    EPISUITE_df['Experimental OH Rate Constants'] = OH_exp
    EPISUITE_df['Predicted Ozone Reaction Rate (cm3/molecule-sec)'] = ozone_predicted
    EPISUITE_df['Experimental Ozone Reaction Rate'] = ozone_exp

    EPISUITE_df.set_index('SMILES', inplace=True)

    return EPISUITE_df


def predictions_QSARToolbox(output_QSARToolbox, endpoint):
    QTB_df = pd.DataFrame()

    smiles_list = []
    endpoint_value = []
    endpoint_AD = []
    mask = output_QSARToolbox['Endpoint']==endpoint
    df_endpoint = output_QSARToolbox[mask][['SMILES', 'Value.MeanValue', 'Domain']]
    df_endpoint['SMILES'] = df_endpoint['SMILES'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)))
    df_endpoint['Response value'] = df_endpoint['Value.MeanValue'].map({'Positive':1, 'Negative':0})
    df_endpoint = df_endpoint.groupby('SMILES').aggregate({'Value.MeanValue':lambda x:igualdad_lista(x), 'Domain': lambda x:igualdad_lista(x)})
    df_endpoint.reset_index(inplace=True)
    for row in df_endpoint.iterrows():
        smiles_list.append(row[1]['SMILES'])
        endpoint_value.append(row[1]['Value.MeanValue'])
        endpoint_AD.append(row[1]['Domain'].split('(')[0].strip())

    QTB_df['SMILES'] = smiles_list
    QTB_df[f'{endpoint} value'] = endpoint_value
    QTB_df[f'{endpoint} AD'] = endpoint_AD

    QTB_df.set_index('SMILES', inplace=True)

    return QTB_df
    

def predictions_ADMETlab(output_ADMETlab):
    ADMETlab_df = pd.DataFrame()

    smiles_list = []
    oral_acute = []
    respiratory = []

    for row in output_ADMETlab.iterrows():
        smiles_list.append(Chem.MolToSmiles(Chem.MolFromSmiles(row[1]['smiles'])))
        if row[1]['ROA'] < 0.5:
            oral_acute.append(0)
        else:
            oral_acute.append(1)

        if row[1]['Respiratory'] < 0.5:
            respiratory.append(0)
        else:
            respiratory.append(1)

    ADMETlab_df['SMILES'] = smiles_list
    ADMETlab_df['Acute oral toxicity'] = oral_acute
    ADMETlab_df['Inhalation toxicity'] = respiratory

    ADMETlab_df.set_index('SMILES', inplace=True)

    return ADMETlab_df

def search_predictions_batch():
    QSAR_df = pd.DataFrame()

    for source in QSAR_MODELS:
        if source == "ProtoPRED":
            output_ProtoTOX = pd.read_excel(OUTPUTS_PATH + ProtoTOX_file, sheet_name='Summary')
            ProtoTOX_df = predictions_ProtoPRED(output_ProtoTOX)
            QSAR_df = QSAR_df.join(ProtoTOX_df, how='outer') 

            output_ProtoECO = pd.read_excel(OUTPUTS_PATH + ProtoECO_file, sheet_name='Summary')
            ProtoECO_df = predictions_ProtoPRED(output_ProtoECO)
            QSAR_df = QSAR_df.join(ProtoECO_df, how='outer')

            output_ProtoPHYSCHEM = pd.read_excel(OUTPUTS_PATH + ProtoPHYSCHEM_file, sheet_name='Summary')
            ProtoPHYSCHEM_df = predictions_ProtoPRED(output_ProtoPHYSCHEM)
            QSAR_df = QSAR_df.join(ProtoPHYSCHEM_df, how='outer')

            output_ProtoADME = pd.read_excel(OUTPUTS_PATH + ProtoADME_file, sheet_name=None)
            prop = ProtoADME_dict(output_ProtoADME)
            ProtoADME_df = predictions_ProtoADME(output_ProtoADME, prop)
            QSAR_df = QSAR_df.join(ProtoADME_df, how='outer')
            
        elif source=='ProtoAquaTox':
            output_ProtoAquaTox = pd.read_excel(OUTPUTS_PATH + ProtoAquaTox_file, engine='openpyxl')

            for endpoint in QSAR_MODELS['ProtoAquaTox']:
                if endpoint=='Acute toxicity on invertebrates':
                    AcToxInv = predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint)
                    QSAR_df = QSAR_df.join(AcToxInv, how='outer')

                elif endpoint=='Acute growth inhibition on algae':
                    AcGrowthAlg = predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint)
                    QSAR_df = QSAR_df.join(AcGrowthAlg, how='outer')

                elif endpoint=='Long-term toxicity on fish':
                    ChronicToxFish = predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint)
                    QSAR_df = QSAR_df.join(ChronicToxFish, how='outer')

                elif endpoint=='Long-term toxicity on invertebrates':
                    ChronicToxInv = predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint)
                    QSAR_df = QSAR_df.join(ChronicToxInv, how='outer')

                elif endpoint=='Chronic growth inhibition on algae':
                    ChronicGrowthAlg = predictions_ProtoAquaTox(output_ProtoAquaTox, endpoint)
                    QSAR_df = QSAR_df.join(ChronicGrowthAlg, how='outer')


        elif source=="Vega":
            output_Vega = pd.read_csv(OUTPUTS_PATH + VEGA_file, sep='\t', skiprows=6, encoding='latin')
            VEGA_df = predictions_Vega(output_Vega)
            
            QSAR_df = QSAR_df.join(VEGA_df, how='outer')

        elif source=='TEST':
            for endpoint in QSAR_MODELS['TEST']:
                TEST_file = OUTPUTS_PATH + endpoint + '_Consensus.csv'
                output_TEST = pd.read_csv(TEST_file, sep=',', encoding='latin')
                TEST_df = predictions_TEST(output_TEST, endpoint)

                QSAR_df = QSAR_df.join(TEST_df, how='outer')

        elif source=='KATE':
            #No me acaba de gustar elmodelo y el output que da. Yo lo quitaría de la lista.
            pass

        elif source=='MLTOX':
            output_MLTOX = pd.read_csv(OUTPUTS_PATH+MLTOX_file, sep=';')
            MLTOX_df = predictions_MLTOX(output_MLTOX)

            QSAR_df = QSAR_df.join(MLTOX_df, how='outer')

        elif source=='EPISUITE':
            
            output_EPISUITE = pd.read_csv(OUTPUTS_PATH+EPISUITE_file, engine='python', sep='\s\s\s*', skiprows=2, names=['OH predicted', 'OH experimental', 'Ozone predicted', 'Ozone experimental', 'SMILES+Id'])
            EPISUITE_df = predictions_EPISUITE(output_EPISUITE)

            QSAR_df = QSAR_df.join(EPISUITE_df, how='outer')

        elif source=='QSARToolbox':
            output_QSARToolbox = pd.read_csv(OUTPUTS_PATH+QTB_file, sep='\t', encoding='utf-16')

            for endpoint in QSAR_MODELS['QSARToolbox']:
                QTB_df = predictions_QSARToolbox(output_QSARToolbox, endpoint)
                QSAR_df = QSAR_df.join(QTB_df, how='outer')

        elif source=='ADMETlab':
            output_ADMETlab = pd.read_csv(OUTPUTS_PATH+ADMETlab_file)
            ADMETlab_df = predictions_ADMETlab(output_ADMETlab)

            QSAR_df = QSAR_df.join(ADMETlab_df, how='outer')
            

        elif source=='STopTox':
            #Mirar con Salva si existía una API
            print("WARNING: No output treatment possible because the platform doesn't allow to predict a batch of compounds")

        else:
            print("WARNING: No output treatment script has been created for models from ", source)
        
    # VEGA_df.to_csv(OUTPUTS_PATH+'VEGA_df_test.csv', sep=';')
    # TEST_df.to_csv(OUTPUTS_PATH+'TEST_df_test.csv', sep=';')
    # MLTOX_df.to_csv(OUTPUTS_PATH+'MLTOX_df_test.csv', sep=';')
    # EPISUITE_df.to_csv(OUTPUTS_PATH+'EPISUITE_df_test.csv', sep=';')
    # for endpoint in QSAR_MODELS['QSARToolbox']:
    #     no_space = endpoint.replace(' ','')
    #     exec(f"{no_space}_df.to_csv(OUTPUTS_PATH+'{no_space}_df_test.csv', sep=';')")
    # ADMETlab_df.to_csv(OUTPUTS_PATH+'ADMETlab_df_test.csv', sep=';')

    QSAR_df.to_csv(OUTPUTS_PATH+'QSAR_df_test.csv', sep=';')


if __name__ == '__main__':
    main()

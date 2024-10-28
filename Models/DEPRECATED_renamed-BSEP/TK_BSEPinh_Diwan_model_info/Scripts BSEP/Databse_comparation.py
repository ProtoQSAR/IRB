#%%
import pandas as pd
from rdkit import Chem

PATH = 'C:/Users/Enrique/Desktop/'

def sanitize(df):
    san_smiles = []
    incorrect_smiles = []

    # Santize smiles
    smiles = list(df['SMILES'])
    for sm in smiles:
        try:
            sanitized_smi = Chem.MolToSmiles(Chem.MolFromSmiles(sm))
            san_smiles.append(sanitized_smi)
        except:
            san_smiles.append('-')
            incorrect_smiles.append(sm)
    
    return san_smiles


def classification_metrics(y_exp, y_pred, y_exp_tag, y_pred_tag, positive_tag, negative_tag):
    
    print(f"{(y_exp == 1).sum()} {positive_tag} and {(y_exp == 0).sum()} {negative_tag} in {y_exp_tag}\n")
    print(f"{(y_pred == 1).sum()} {positive_tag} and {(y_pred == 0).sum()} {negative_tag} in {y_pred_tag}\n")
    
    FP=0
    FN=0
    TP=0
    TN=0

    for i in range(len(y_exp)):
        if y_exp[i] == 0 and y_pred[i] == 0:
            TN=TN+1
        elif y_exp[i] == 1 and y_pred[i] == 0:
            FN=FN+1
        elif y_exp[i] == 1 and y_pred[i] == 1:
            TP=TP+1
        elif y_exp[i] == 0 and y_pred[i] == 1:
            FP=FP+1

    # Keep this print to chek your values
    print(f"\n\n{y_exp_tag} vs {y_pred_tag}")
    print(TP,TN,FP,FN)
    
    print('\n\t\t\t\tPredcited')
    print(f'\t\t\t\t{negative_tag}\t{positive_tag}')
    print(f'Actual\t{negative_tag}\t'+str(TN)+'\t'+str(FP))
    print(f'\t\t{positive_tag}\t'+str(FN)+'\t\t'+str(TP))
    
    accuracy = (TP + TN) / (TP + TN + FP + FN)
    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    precision = TP / (TP + FP)
    f1 = (2 * sensitivity * precision) / (sensitivity + precision)
    
    accuracy,sensitivity,specificity,precision,f1 = round(accuracy,2),round(sensitivity,2),round(specificity,2),round(precision,2),round(f1,2)
    
    # Keep this print to chek your values
    print('\nAccuracy:',accuracy)
    print('Sensitivity:',sensitivity)
    print('Specificity:',specificity)
    print('Precision:',precision)
    print('f1:',f1)
    
    
    return [y_exp_tag, y_pred_tag, len(y_exp), accuracy,sensitivity,specificity,precision,f1, TN,FP,FN,TP]

def compare_compounds(df1, df2):
    """
    Compara los compuestos en dos DataFrames de pandas, e imprime los compuestos que son similares entre ambos,
    verificando además si el valor de la variable 'y' es el mismo para estos compuestos similares.

    Parámetros:
    - df1: Primer DataFrame de pandas que contiene columnas con compuestos (SMILES) y sus valores de 'y'.
    - df2: Segundo DataFrame de pandas que contiene columnas con compuestos (SMILES) y sus valores de 'y'.
    """
    compuestos_df1 = df1.set_index('Sanitized_SMILES')['y'].to_dict()
    compuestos_df2 = df2.set_index('Sanitized_SMILES')['y'].to_dict()

    # Encontrar la intersección de ambos conjuntos para identificar compuestos similares
    compuestos_similares = set(compuestos_df1.keys()).intersection(compuestos_df2.keys())

    suma = 0
    coincidencias = []
    
    if compuestos_similares:
        print("Compuestos similares encontrados:", len(compuestos_similares))
        for compuesto in compuestos_similares:
            valor_y_df1 = compuestos_df1[compuesto]
            valor_y_df2 = compuestos_df2[compuesto]
            coincidencias.append([compuesto, valor_y_df1, valor_y_df2])
            if valor_y_df1 != valor_y_df2:
                suma += 1
                print(f"Compuesto: {compuesto}, Valor 'y' en Validacion: {valor_y_df1}, Valor 'y' en Prediccion: {valor_y_df2}")
    else:
        print("No se encontraron compuestos similares.")

    print("-" * 79)
    print(suma, "compuestos no coinciden en sus valores 'y'")

    # Crear un nuevo DataFrame con las coincidencias
    coincidencias_df = pd.DataFrame(coincidencias, columns=['SMILES', 'y_exp', 'y_pred'])

    return coincidencias_df

def main():

    # Importa las bases de datos
    df1 = pd.read_csv(PATH+'smiles BSEP.txt', sep='\t')
    df2 = pd.read_csv(PATH+'BSEP_experimental.txt', sep='\t')

    smile_list1 = sanitize(df1)
    smile_list2 = sanitize(df2)

    df1['Sanitized_SMILES'] = smile_list1
    df2['Sanitized_SMILES'] = smile_list2

    # Selecciona y renombra las columnas
    df1 = df1[['Sanitized_SMILES', 'y']]
    df2 = df2[['Sanitized_SMILES', 'y']]

    # Procesamiento específico para la base de datos 2
    df1['y_processed'] = df1['y'].apply(lambda x: x.replace('<', '').replace('>', ''))
    df1['y'] = df1['y_processed'].apply(lambda x: 1 if float(x) < 25 else 0 if float(x) > 100 else None)
    df1 = df1.dropna()

    # Comparar compuestos
    df_new = compare_compounds(df1, df2)
    
    # Obtener métricas    
    y_exp = df_new['y_exp']
    y_pred = df_new['y_pred']
    y_exp_tag = 'Experimental'
    y_pred_tag = 'Precicted'
    positive_tag = 'Inhibitor'
    negative_tag = 'Non-Inhibitor'

    classification_metrics(y_exp, y_pred, y_exp_tag, y_pred_tag, positive_tag, negative_tag)

if __name__ == "__main__":
    main()
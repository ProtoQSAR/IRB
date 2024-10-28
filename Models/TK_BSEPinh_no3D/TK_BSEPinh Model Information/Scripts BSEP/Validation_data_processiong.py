#%%
import pandas as pd


PATH = 'C:/Users/Enrique/Desktop/'
validation_df = pd.read_csv(PATH+'smiles BSEP.txt', sep='\t')

# Procesar la columna del endpoint
validation_df['y'] = validation_df['y'].apply(lambda x: x.replace('<','').replace('>',''))
validation_df['y'] = validation_df['y'].apply(lambda x: 1 if float(x) < 25 else 0 if float(x) > 100 else None)

validation_df = validation_df.dropna()

validation_df.to_csv('validation.csv', sep=';', index=False)


#%%
import pandas as pd
validation_df = pd.read_excel('ELlobet_ProtoADME_20240627_103108.xlsx', sheet_name='Bsep nod3d')
validation_df['Chemical name'] = validation_df['Chemical name'].apply(lambda x: 1 if float(x) == 1.0 else 0)
validation_df['Experimental value*'] = validation_df['Experimental value*'].apply(lambda x: 1 if x == 'Inhibitor' else 0 if x == 'Non-inhibitor' else None)
validation_df['Experimental value*'] = validation_df['Experimental value*'].apply(lambda x: 1 if float(x) == 1.0 else 0)

validation_df = validation_df[['Chemical name', 'Experimental value*']]

validation_df = validation_df.dropna()

validation_df.to_csv('Comparativa_databases_BSEP.csv', sep=';', index=False)

#%%
import pandas as pd
validation_df = pd.read_csv('AID_1449628_datatable_1uM.csv') 
# validation_df[['PUBCHEM_EXT_DATASOURCE_SMILES', 'PubChem Standard Value']]

validation_df = validation_df.loc[5:, :]
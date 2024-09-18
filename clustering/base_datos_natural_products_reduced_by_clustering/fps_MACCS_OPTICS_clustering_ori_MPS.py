
##
#%%############################# CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i)
##############################################################################

#%% imports
import pandas as pd
import os
import numpy as np

from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn.cluster import OPTICS

import pickle

# Script to calculate MACCSkeys fingerprints of molecules, 
#save them by chunks in .sav files, as objects to save space
#and then use the optics algorithm to cluster them and 
#perform a subselection of molecules

#%% functions

def fps_maccs_daylight(lista):
    # tanim_topol = []
    tanim_maccs = []
    for i, molec in enumerate(lista):
        # print(f"Molecule number {i}: {molec}")
        ref = Chem.MolFromSmiles(molec)
        try:
            # fp = Chem.RDKFingerprint(ref)
            fp_macc = MACCSkeys.GenMACCSKeys(ref)
            
            # tanim_topol.append(fp)
            tanim_maccs.append(fp_macc)
        except:
            print(f"{molec} gives an error")
        
    # return tanim_topol,tanim_maccs
    return tanim_maccs

#fingerprints by chunks save to parquet (less memory)    
def fps_to_pickle(df_in, splits):

    for i,df in enumerate(np.array_split(df_in, splits)):
            
        chunk_name = f"chunk{str(i)}"
        print(chunk_name, df.shape)
        smis = list(df.SMILES)
        
        print("calculating fingerprints...")
        fps_maccs = fps_maccs_daylight(smis)
        
        df["MACCS"] = fps_maccs
        outfile = chunk_name+'_maccs_fps.sav'
        pickle.dump(df, open(outfile, "wb"))
        
        print(f"Fps saved in file {outfile}") 
        
def fps_array(pckld_df):
    fps = pckld_df.MACCS
    fps_arr = np.array([list(fp) for fp in fps])
    return fps_arr 
		
		
#%%read files
####IMPORTANT: Delete the generate files before running the script again
data = "tox21-er-bla-agonist-p2-string_not_included_NURA"
path = 'C:/Users/Martina/Documents/Endocrine_disruptors/Modelos_estrogen_receptor/Model_AGO_ERA_2class/Validation/tox21/'
os.chdir(path)
study_frame = pd.read_csv (path+f"{data}.csv", sep=';', encoding='latin-1')

		
#%%convert float to string
####IMPORTANT: Delete the generate files before running the script again
study_frame['SMILES'] = study_frame['SMILES'].astype(str) 
study_frame.to_csv("tox21-er-bla-agonist-p2-string-preprocessed.csv", sep=";",index=False)   

		
#%%Create a subset with only the SMILES and the y columns
study_frame_sub = study_frame.loc[:,["SMILES","y"]]

		
#%%Create a subset with only y= 0 values
study_frame_sub = study_frame_sub[study_frame_sub["y"]==0]
print(study_frame_sub.columns)
		
#%%Divide the dataset in parts (you can change the number of parts below)

fps_to_pickle(study_frame_sub, 100) 

#%% read savs and cluster with optics
files = os.listdir()
sav_files = [file  for file in files if file.endswith("sav")]

df_all = pd.DataFrame()
for i,file in enumerate(sav_files):
    if i == 0:
        df_all = pickle.load(open(file, "rb"))
    else:
        df_sav = pickle.load(open(file, "rb"))
        df_all = pd.concat([df_all,df_sav])

df_all.to_csv(f"{str(df_all.shape[0])}molecs_wMACCSfps.csv", sep=";",index=False)
        
arr_fps = fps_array(df_all)     
    
#Clustering with OPTICS
 # n_jobs t parallelize     
clustering = OPTICS(metric="jaccard", n_jobs=-12).fit(arr_fps)  
labels = clustering.labels_

#%% save df of all molecules with cluster tags
    
df_all["clusters"] = labels
    
df_all.to_csv(f"{str(df_all.shape[0])}molecs_wMACCSfps_clusters.csv", sep=";",index=False)    

#%% create dictionary of clusters    
cluster_dic = {}
for i in df_all.clusters.unique():
    cluster_name = f"cluster_{str(i)}"
    cluster_dic[cluster_name] = df_all[df_all["clusters"] == i]

outfile = f'{str(len(cluster_dic))}_clusters_dictionary.sav'
pickle.dump(cluster_dic, open(outfile, "wb"))

#%% subset of each cluster (change frac to chose the number of molecules)
subsets_dic = {}
for key, df in cluster_dic.items():
    subsets_dic[key] = df.sample(frac=0.15)

df_sub_all = pd.concat(subsets_dic.values(), ignore_index=True) 
df_sub_all.to_csv(f"{str(df_sub_all.shape[0])}molecs_subset_clusters.csv", sep=";",index=False)   


#%% generate final dataset with selected compounds

   # Filter values y =  0 of initial dataset
# study_frame_0 = study_frame[df['Cluster'].isin([1,2, 3]) & (df['y'] != 1)]
study_frame_0 = study_frame[(study_frame['y'] == 1)]

# Merge selected smiles with initial value
selected_smiles = df_sub_all['SMILES']
study_frame_1 = study_frame['SMILES'].isin(selected_smiles)
study_frame_1 = study_frame[study_frame['SMILES'].isin(selected_smiles)]
# Concat dataframes with and without weakly active
df_final = pd.concat([study_frame_0, study_frame_1])
   
df_final.to_csv(path+f"{data}_0selected.csv", sep=";",index=False)   

# %%

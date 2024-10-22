"""
Description: A Script to perform KMeans clustering on a dataset of molecules using MACCS fingerprints.
It reads a CSV file with SMILES strings, calculates MACCS fingerprints, and saves the fingerprints in a pickle file.
Then, it reads the pickle files, performs KMeans clustering for different numbers of clusters, and plots the inertia values
to help identify the optimal number of clusters using the elbow method. Additionally, it prints the silhouette score for each clustering attempt. 
This is based on different functions developed by Rita Ortega (PQS) for OPTICS and was adapted for Kmeans by Laureano E. Carpio (Moldrug).

Author(s): Rita Ortega (PQS), Laureano E. Carpio (Moldrug)

Date: 2024-10-16

"""
#%% imports and functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from rdkit.Chem import MACCSkeys
import pickle
import os
from rdkit import Chem
from sklearn.metrics import silhouette_score

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



def optimal_kmeans_elbow(arr_fp, first_number = 2, max_clusters=10, step=1):
    """
    Perform KMeans clustering on an array for different cluster numbers,
    and plot the inertia to help identify the optimal number of clusters using the elbow method.
    Additionally, prints the silhouette score for each clustering attempt.

    Parameters:
    - arr_fp: Array or DataFrame with the data for clustering
    - max_clusters: Maximum number of clusters to test (default is 10)

    Returns:
    - A plot of inertia values for different cluster numbers
    - The optimal number of clusters (integer) based on the elbow method
    """
    # Store the inertia values and silhouette scores for each k
    inertia_values = []
    silhouette_scores = []

    # Test KMeans for each number of clusters
    for k in range(first_number, max_clusters + 1, step):  # Start from 2 clusters, as silhouette score is undefined for 1

        print(f" [+] PROCESS WITH {k} CLUSTER(S)")
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(arr_fp)
        inertia_values.append(kmeans.inertia_)

        # Calculate silhouette score
        silhouette_avg = silhouette_score(arr_fp, kmeans.labels_)
        silhouette_scores.append(silhouette_avg)
        print(f"    Silhouette Score for {k} clusters: {silhouette_avg}")

    # Plot the inertia values for each number of clusters
    plt.figure(figsize=(10, 6))
    plt.plot(range(first_number, max_clusters + 1,step), inertia_values, marker='o')
    plt.title('Elbow Method for Optimal k')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Inertia')
    plt.show()

    # Find the "elbow" point using the inertia values
    optimal_clusters = np.argmin(np.diff(inertia_values, 2)) + 2

    print(f"Optimal number of clusters based on the elbow method: {optimal_clusters}")

    #Same plot but for silhouette scores
    plt.figure(figsize=(10, 6))
    plt.plot(range(first_number, max_clusters + 1,step), silhouette_scores, marker='o')
    plt.title('Silhouette Score for Optimal k')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Silhouette Score')
    plt.show()


    return optimal_clusters

def fps_array(pckld_df):
    fps = pckld_df.MACCS
    fps_arr = np.array([list(fp) for fp in fps])
    return fps_arr

#%%read files
####IMPORTANT: Delete the generate files before running the script again
data = "IRB_library_merged_all"
path = '/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/kmeans_results/'
os.chdir(path)
study_frame = pd.read_csv (path+f"{data}.csv", sep=';', encoding='latin-1')
#%%convert float to string
####IMPORTANT: Delete the generate files before running the script again
study_frame['SMILES'] = study_frame['SMILES'].astype(str)
study_frame.to_csv(data+"-prepared.csv", sep=";",index=False)


#%%Create a subset with only the SMILES and the y columns
#study_frame_sub = study_frame.loc[:,["SMILES","y"]]
study_frame_sub = study_frame

#%%Create a subset with only y= 0 values
#study_frame_sub = study_frame_sub[study_frame_sub["y"]==0]

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

# %%

optim_clu = optimal_kmeans_elbow(arr_fps, first_number=550, max_clusters=1050, step = 50)

# %%

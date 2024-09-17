# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:10:04 2022

@author: Rita
"""

import os 
import pandas as pd
import pickle
import numpy as np
from sklearn.cluster import OPTICS

# os.chdir("F:/RITA/Practicas/base_datos_Laureano")

# df_sav = pickle.load(open("FINAL_429_maccs_fps.sav", "rb"))

os.chdir(r"C:/Users/proto/Desktop/Rita/base_datos_Laureano/1intento_clustering")
#%% build array of fingerprints

def fps_array(pckld_df):
    fps = pckld_df.MACCS
    fps_arr = np.array([list(fp) for fp in fps])
    return fps_arr 

# array = fps_array(df_sav)
# xktin = array[:100,:]

#%% clustering

clustering = OPTICS(metric="jaccard").fit(xktin)
print(len(clustering.labels_))

labels = clustering.labels_

df_mini = df_sav.iloc[:100,:]
df_mini["clusters"] = labels

# n_clusters = len(df_mini.clusters.unique())

cluster_dic = {}
for i in df_mini.clusters.unique():
    cluster_name = f"cluster_{str(i)}"
    cluster_dic[cluster_name] = df_mini[df_mini["clusters"] == i]

subsets_dic = {}
for key, df in cluster_dic.items():
    subsets_dic[key] = df.sample(frac=0.1)

df_sub_all = pd.concat(subsets_dic.values(), ignore_index=True)

#%% functions to apply to all the dataset

files = os.listdir()
sav_files = [file  for file in files if file.endswith("sav")]

df_all = pd.DataFrame()
for i,file in enumerate(sav_files):
    if i == 0:
        df_all = pickle.load(open(file, "rb"))
    else:
        df_sav = pickle.load(open(file, "rb"))
        df_all = pd.concat([df_all,df_sav])
        
arr_fps = fps_array(df_all)  
clustering = OPTICS(metric="jaccard").fit(arr_fps)  
labels = clustering.labels_

df_all["clusters"] = labels

cluster_dic = {}
for i in df_all.clusters.unique():
    cluster_name = f"cluster_{str(i)}"
    cluster_dic[cluster_name] = df_all[df_all["clusters"] == i]

subsets_dic = {}
for key, df in cluster_dic.items():
    subsets_dic[key] = df.sample(frac=0.1)

df_sub_all = pd.concat(subsets_dic.values(), ignore_index=True)
    
outfile = f'{str(df_sub_all.shape[0])}molecs_subset_clusters.sav'
pickle.dump(df_sub_all, open(outfile, "wb"))
    
df_sub_all.to_csv(f"{str(df_sub_all.shape[0])}molecs_subset_clusters.csv", sep=";",index=False)   

#%% try another clustering and save DICT of clusters

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

###
    
df_all["clusters"] = labels
    
df_all.to_csv(f"{str(df_all.shape[0])}molecs_wMACCSfps_clusters.csv", sep=";",index=False)    
    
cluster_dic = {}
for i in df_all.clusters.unique():
    cluster_name = f"cluster_{str(i)}"
    cluster_dic[cluster_name] = df_all[df_all["clusters"] == i]

outfile = f'{str(len(cluster_dic))}_clusters_dictionary.sav'
pickle.dump(cluster_dic, open(outfile, "wb"))

subsets_dic = {}
for key, df in cluster_dic.items():
    subsets_dic[key] = df.sample(frac=0.1)

df_sub_all = pd.concat(subsets_dic.values(), ignore_index=True) 
df_sub_all.to_csv(f"{str(df_sub_all.shape[0])}molecs_subset_clusters.csv", sep=";",index=False)   
   

#%% imports and functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
# from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
import pickle
import os
from rdkit import Chem
from sklearn.metrics import silhouette_score
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
from sklearn.cluster import OPTICS
from sklearn.metrics import silhouette_score
from sklearn.cluster import DBSCAN
from sklearn.cluster import HDBSCAN

def fps_array(pckld_df):
    fps = pckld_df.Morgan
    fps_arr = np.array([list(fp) for fp in fps])
    return fps_arr

def fps_maccs_daylight(lista):
    # tanim_topol = []
    tanim_morgan = []
    for i, molec in enumerate(lista):
        # print(f"Molecule number {i}: {molec}")
        ref = Chem.MolFromSmiles(molec)
        try:
            # fp = Chem.RDKFingerprint(ref)
            # fp_macc = MACCSkeys.GenMACCSKeys(ref)
            fp_morgan_prev = AllChem.GetMorganFingerprintAsBitVect(ref, useChirality=True, radius=3, nBits = 1024)
            fp_morgan = np.array(fp_morgan_prev)

            # tanim_topol.append(fp)
            tanim_morgan.append(fp_morgan)
        except:
            print(f"{molec} gives an error")

    # return tanim_topol,tanim_maccs
    return tanim_morgan

def fps_to_pickle(df_in, splits):

    for i,df in enumerate(np.array_split(df_in, splits)):

        chunk_name = f"chunk{str(i)}"
        print(chunk_name, df.shape)
        smis = list(df.SMILES)

        print("calculating fingerprints...")
        fp_morgan = fps_maccs_daylight(smis)

        df["Morgan"] = fp_morgan
        outfile = chunk_name+'_morgan_fps.sav'
        pickle.dump(df, open(outfile, "wb"))

        print(f"Fps saved in file {outfile}") 

def pca_analysis(fingerprints, variance_threshold=0.80):
    # Perform PCA
    pca = PCA()
    pca.fit(fingerprints)
    
    # Calculate the cumulative explained variance
    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance_ratio)
    
    # Plot explained variance
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker='o', linestyle='--')
    plt.axhline(y=variance_threshold, color='r', linestyle='-')
    plt.title('Explained Variance by Principal Components')
    plt.xlabel('Number of Principal Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.grid(True)
    plt.show()
    
    # Find the number of components that explain the desired variance
    num_components = np.argmax(cumulative_variance >= variance_threshold) + 1
    print(f'Number of components that explain {variance_threshold*100}% of the variance: {num_components}')
    
    # Perform PCA with the selected number of components
    pca = PCA(n_components=num_components)
    principal_components = pca.fit_transform(fingerprints)
    
    return principal_components

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
        # silhouette_avg = silhouette_score(arr_fp, kmeans.labels_)
        # silhouette_scores.append(silhouette_avg)
        # print(f"    Silhouette Score for {k} clusters: {silhouette_avg}")

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
    # plt.figure(figsize=(10, 6))
    # plt.plot(range(first_number, max_clusters + 1,step), silhouette_scores, marker='o')
    # plt.title('Silhouette Score for Optimal k')
    # plt.xlabel('Number of Clusters (k)')
    # plt.ylabel('Silhouette Score')
    # plt.show()
#%%read files
####IMPORTANT: Delete the generate files before running the script again
data = "IRB_library_merged_all"
path = '/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/'
os.chdir(path)
study_frame = pd.read_csv (path+f"{data}.csv", sep=';', encoding='latin-1')
#%%convert float to string
####IMPORTANT: Delete the generate files before running the script again
study_frame['SMILES'] = study_frame['SMILES'].astype(str)
study_frame.to_csv(data+"-prepared.csv", sep=";",index=False)


#%%Create a subset with only the SMILES and the y columns
#study_frame_sub = study_frame.loc[:,["SMILES","y"]]
study_frame_sub = study_frame
del study_frame_sub["y"]
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

df_all.to_csv(f"{str(df_all.shape[0])}molecs_MORGAN-fps.csv", sep=";",index=False)

arr_fps = fps_array(df_all)
# %%
# principal_components = pca_analysis(arr_fps, variance_threshold=0.80)
# %%
# optim_clu = optimal_kmeans_elbow(principal_components, first_number=20000, max_clusters=50000, step = 5000)
print("[++] Comienza clustering con 40k grupos")
kmeans = KMeans(n_clusters=40000, random_state=42)
kmeans.fit(arr_fps)
#Save the model
pickle.dump(kmeans, open("kmeans_model_40k.sav", "wb"))
print("[++] Acaba clustering con 40k grupos")


print("[+++++] Comienza el cálculo de métricas con clustering con 40k grupos")
# print the silhouette score
silhouette_avg = silhouette_score(arr_fps, kmeans.labels_)

#print the label number and the number of entries in each cluster
# print("Number of clusters: ", len(set(kmeans.labels_)))
print("Silhouette Score: ", silhouette_avg)
print("Number of entries in each cluster: ")
print(pd.Series(kmeans.labels_).value_counts())

print("[+++++] Acaba el cálculo de métricas con clustering con 40k grupos")
#%%
print("[++] Comienza clustering con 45k grupos")
kmeans = KMeans(n_clusters=45000, random_state=42)
kmeans.fit(arr_fps)
#Save the model
pickle.dump(kmeans, open("kmeans_model_45k.sav", "wb"))

print("[++] Acaba clustering con 45k grupos")


print("[+++++] Comienza el cálculo de métricas con clustering con 45k grupos")

# print the silhouette score
silhouette_avg = silhouette_score(arr_fps, kmeans.labels_)

#print the label number and the number of entries in each cluster
# print("Number of clusters: ", len(set(kmeans.labels_)))
print("Silhouette Score: ", silhouette_avg)
print("Number of entries in each cluster: ")
print(pd.Series(kmeans.labels_).value_counts())

print("[+++++] Acaba el cálculo de métricas con clustering con 45k grupos")

#%%
print("[++] Comienza clustering con 50k grupos")
kmeans = KMeans(n_clusters=50000, random_state=42)
kmeans.fit(arr_fps)
#Save the model
pickle.dump(kmeans, open("kmeans_model_50k.sav", "wb"))
print("[++] Acaba clustering con 50k grupos")

print("[+++++] Comienza el cálculo de métricas con clustering con 50k grupos")

# print the silhouette score
silhouette_avg = silhouette_score(arr_fps, kmeans.labels_)

#print the label number and the number of entries in each cluster
# print("Number of clusters: ", len(set(kmeans.labels_)))
print("Silhouette Score: ", silhouette_avg)
print("Number of entries in each cluster: ")
print(pd.Series(kmeans.labels_).value_counts())

print("[+++++] Acaba el cálculo de métricas con clustering con 50k grupos")
#%% Count the number of clusters with less than 3 entries
aa = pd.Series(kmeans.labels_).value_counts()
len(aa[aa>=3])
morethan3 = aa[aa>=3]
sum(morethan3)
# %% Add the cluster labels to the dataframe df_all in a column called Cluster_first_round
kmeans_40k = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/kmeans_model_40k.sav', 'rb'))

labels_40k = kmeans_40k.labels_

df_all["Cluster"] = labels_40k
df_all = df_all.drop(columns=["Morgan"])
df_all.to_csv("IRB_library_merged_clustering_50k.csv", sep=";",index=False)


#%% Add a new column innicating the center of each of the clusters

import numpy as np
import pandas as pd

# Crear una lista inicializada con False para almacenar si una entrada es el centroide
is_centroid = np.zeros(len(df_all), dtype=bool)

# Iterar sobre cada clúster
for cluster_id in range(kmeans.n_clusters):
    # Obtener los índices de las entradas que pertenecen al clúster actual
    cluster_indices = np.where(df_all['Cluster'] == cluster_id)[0]
    cluster_points = arr_fps[cluster_indices]  # Seleccionar los puntos del clúster

    # Si no hay puntos en el clúster, pasar al siguiente
    if len(cluster_points) == 0:
        continue

    # Obtener el centroide del clúster
    centroid = kmeans.cluster_centers_[cluster_id]

    # Calcular la distancia de cada punto al centroide
    distances = np.linalg.norm(cluster_points - centroid, axis=1)

    # Identificar el índice relativo del punto más cercano al centroide
    closest_index_relative = np.argmin(distances)

    # Convertir el índice relativo al índice absoluto
    closest_index_absolute = cluster_indices[closest_index_relative]

    # Marcar únicamente esta entrada como el centroide
    is_centroid[closest_index_absolute] = True

# Agregar la columna al DataFrame
df_all['is_centroid'] = is_centroid


#%% sort in ascending order based on the cluster column
df_all2 = df_all.sort_values(by="Cluster")

#%% save the fileç
df_all2.to_csv("IRB_library_merged_clustering40K_sorted.csv", sep=";",index=False)
# #%%
df_all2.to_excel("IRB_library_merged_clustering40K_sorted.xlsx",index=False)

# #%%
df_all.to_csv("IRB_library_merged_clustering40K_not-sorted.csv", sep=";",index=False)
# #%%
df_all.to_excel("IRB_library_merged_clustering40K_not-sorted.xlsx",index=False)

# %%count the number of clusters with less than 5 entries
len(aa[aa<5])

# %%count the number of compounds that pertain to clusterswith less of 3 entries
lessthan3 = aa[aa<3]
# # sum all the values
sum(lessthan3)

# #%%
morethan3 = aa[aa>=3]
sum(morethan3)


# %%
# import numpy as np
# from scipy.spatial.distance import cdist

# # # Obtener los centroides del modelo
# centroids = kmeans.cluster_centers_

# # Calcular la matriz de distancias entre los centroides
# distance_matrix = cdist(centroids, centroids, metric='euclidean')
# %%
import numpy as np
from scipy.spatial.distance import cdist

# Obtener los centroides del modelo
centroids = kmeans.cluster_centers_

# Tamaño del batch (ajústalo según la memoria)
batch_size = 2000

# Número de centroides
num_centroids = len(centroids)

# Crear una matriz vacía para almacenar las distancias
distance_matrix = np.zeros((num_centroids, num_centroids))

# Iterar por batches
for i in range(0, num_centroids, batch_size):
    for j in range(0, num_centroids, batch_size):
        # Seleccionar submatrices (batches) de centroides
        batch_i = centroids[i:i + batch_size]
        batch_j = centroids[j:j + batch_size]
        
        # Calcular distancias entre los batches
        distances = cdist(batch_i, batch_j, metric='euclidean')
        
        # Guardar las distancias calculadas en la matriz completa
        distance_matrix[i:i + batch_size, j:j + batch_size] = distances

# Convertir a DataFrame (opcional)
centroid_distances_df = pd.DataFrame(distance_matrix, 
                                     index=[f'Cluster_{i}' for i in range(num_centroids)],
                                     columns=[f'Cluster_{i}' for i in range(num_centroids)])

# Mostrar las distancias
centroid_distances_df

# %%

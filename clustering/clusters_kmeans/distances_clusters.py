import numpy as np
from scipy.spatial.distance import cdist
import pickle
import pandas as pd

#Load the kmeans model

kmeans = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/new__cluster_muchos/kmeans_model_40k.sav', 'rb'))

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

# guardar el dataframe como csv
centroid_distances_df.to_csv('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/new__cluster_muchos/centroid_distances.csv', sep = ";", index = False)

# %% Load the csv file and plot it as a heatmap
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Load the csv file
centroid_distances_df = pd.read_csv('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/new__cluster_muchos/centroid_distances.csv', sep = ";")

#%%
# Plot the heatmap
sns.heatmap(centroid_distances_df, cmap='viridis', square=True, cbar_kws={'label': 'Euclidean distance'})
plt.title('Distances between cluster centroids')
plt.show()

# %%

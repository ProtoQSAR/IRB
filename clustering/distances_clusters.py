#%%

import numpy as np
from scipy.spatial.distance import cdist
import pickle
import pandas as pd
from openpyxl import Workbook

# Load the kmeans model
kmeans = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/kmeans_model_40k.sav', 'rb'))

# Obtener los centroides del modelo
centroids = kmeans.cluster_centers_

# Tamaño del batch (ajústalo según la memoria)
batch_size = 4000

# Número de centroides
num_centroids = len(centroids)

# Crear una matriz vacía para almacenar las distancias
distance_matrix = np.zeros((num_centroids, num_centroids))

# Ruta base para guardar los archivos Excel
output_dir = '/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/'
output_base_name = 'centroid_distances_batch'

# Iterar por batches
print("Guardando distancias en archivos Excel por batches...")
for i in range(32000, num_centroids, batch_size):
    for j in range(0, num_centroids, batch_size):
        # Seleccionar submatrices (batches) de centroides
        batch_i = centroids[i:i + batch_size]
        batch_j = centroids[j:j + batch_size]
        
        # Calcular distancias entre los batches
        distances = cdist(batch_i, batch_j, metric='euclidean')
        
        # Guardar las distancias calculadas en la matriz completa
        distance_matrix[i:i + batch_size, j:j + batch_size] = distances

        # Convertir a DataFrame
        batch_df = pd.DataFrame(
            distances,
            index=[f'Cluster_{x}' for x in range(i, i + len(batch_i))],
            columns=[f'Cluster_{x}' for x in range(j, j + len(batch_j))]
        )

        # Nombre del archivo Excel para este batch
        output_file = f"{output_dir}{output_base_name}_i{i}_j{j}.xlsx"

        # Guardar el batch como archivo Excel
        batch_df.to_excel(output_file, index=True)
        print(f"Batch guardado: {output_file}")

# Mensaje final
print("Distancias guardadas en archivos Excel por batches.")

# %%

# %%

###############################################################################
############################### IMPORT LIBRARIES ##############################
###############################################################################

import pandas as pd
import numpy as np
import os

from scipy.spatial import distance
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import Descriptors

###############################################################################
################################ PATH AND FILE ################################
###############################################################################

intro = r'C:\Users\Agata\Documents\GitHub\HYPIEND\07Prediction_models\04NEO\3D_binary\ANT-TRA'
modified_intro = intro.replace('\\', '/')
modified_intro += '/'

PATH = modified_intro
name = "ANT-TRA"
df = pd.read_csv(PATH+name+"-initial_reduction.csv", sep=";")

###############################################################################
############################# OUTLIERS ELIMINATION ############################
###############################################################################

# unique_values = df['y'].unique()

# if set(unique_values) == {0, 1}:
#     colors = ['tab:blue', 'tab:red']
# elif set(unique_values) == {0, 1, 2}:
#     colors = ['tab:blue', 'tab:red', 'tab:green']
# else:
#     print("Column 'y' has no valid values")

# # Calculate Molecular Weight
# df['MolecularWeight'] = df['SMILES'].apply(lambda x: Chem.Descriptors.MolWt(Chem.MolFromSmiles(x)))

# # Calculate quartiles
# q3 = df['MolecularWeight'].quantile(0.75)
# q1 = df['MolecularWeight'].quantile(0.25)
# iqr = q3 - q1

# # Calculate upper and lower limits using interquartile range
# upper_limit = q3 + 3 * iqr
# lower_limit = q1 - 3 * iqr

# # Filter outliers
# filtered_outliers = df[(df['MolecularWeight'] > upper_limit) | (df['MolecularWeight'] < lower_limit)]
# filtered_outliers.to_csv(PATH+f"{name}_outliers.csv", sep=";", index=False)

# # Plot box-whisker plot with outliers
# sns.set(rc={'figure.figsize':(20,15)})
# ax = sns.boxplot(data=df, x="y", y="MolecularWeight", palette=colors)

# # Highlight outliers
# ax.scatter(filtered_outliers['y'], filtered_outliers['MolecularWeight'], color='red', marker='o', s=100, label='Outliers')

# plt.title(f"{name}_highlight_outliers")
# plt.legend()
# plt.show()

# # Remove outliers
# remove_outliers = df[(df['MolecularWeight'] >= lower_limit) & (df['MolecularWeight'] <= upper_limit)]
# remove_outliers.to_csv(PATH+f"{name}_removed_outliers.csv", sep=";", index=False)
# df = remove_outliers

# # Plot box-whisker plot without outliers
# sns.set(rc={'figure.figsize':(20,15)})
# sns.boxplot(data=remove_outliers, x="y", y="MolecularWeight", palette=colors).set_title(f"{name}_removed_outliers")

# plt.show()

##############################################################################
############################ VALUES ELIMINATION ##############################
##############################################################################

# df = pd.read_csv(PATH+name+"_removed_outliers.csv", sep=";")
# df = df.loc[df['y'] != 1]
# # df = df.loc[df['y'] != 2]

# ###############################################################################
# ########################## WEIRD VALUES ELIMINATION ###########################
# ###############################################################################

# value_1 = 'C=CC[N+]12CCC34c5ccccc5N5C=C6C7CC8C9(CC[N+]8(CC=C)CC7=CCO)c7ccccc7N(C=C(C(CC31)C(=CCO)C2)C54)C69'
# # value_2 = 'CNC(CC(C)C)C(=O)NC1C(=O)NC(CC(N)=O)C(=O)NC2C(=O)NC3C(=O)NC(C(=O)NC(C(=O)O)c4cc(O)cc(O)c4-c4cc3ccc4O)C(O)c3ccc(c(Cl)c3)Oc3cc2cc(c3OC2OC(CO)C(O)C(O)C2OC2CC(C)(N)C(O)C(C)O2)Oc2ccc(cc2Cl)C1O'
# # value_3 = 'OCC1OC2OC3C(CO)OC(OC4C(CO)OC(OC5C(CO)OC(OC6C(CO)OC(OC7C(CO)OC(OC8C(CO)OC(OC9C(CO)OC(OC1C(O)C2O)C(O)C9O)C(O)C8O)C(O)C7O)C(O)C6O)C(O)C5O)C(O)C4O)C(O)C3O'
# # value_4 = "OCC1OC2OC3C(CO)OC(OC4C(CO)OC(OC5C(CO)OC(OC6C(CO)OC(OC7C(CO)OC(OC8C(CO)OC(OC1C(O)C2O)C(O)C8O)C(O)C7O)C(O)C6O)C(O)C5O)C(O)C4O)C(O)C3O"
# # value_5 = "Cc1ccc(C(=O)Nc2ccc(S(=O)(=O)O)c3cc(S(=O)(=O)O)cc(S(=O)(=O)O)c23)cc1NC(=O)c1cccc(NC(=O)Nc2cccc(C(=O)Nc3cc(C(=O)Nc4ccc(S(=O)(=O)O)c5cc(S(=O)(=O)O)cc(S(=O)(=O)O)c45)ccc3C)c2)c1"

# df = df[df['SMILES'] != value_1].reset_index(drop=True)
# # df = df[df['SMILES'] != value_2].reset_index(drop=True)
# # df = df[df['SMILES'] != value_3].reset_index(drop=True)
# # df = df[df['SMILES'] != value_4].reset_index(drop=True)
# # df = df[df['SMILES'] != value_5].reset_index(drop=True)

##############################################################################
######################### ELBOW METHOD FOR K-MEANS ###########################
##############################################################################

# print("Applying the Elbow Method...")

# # Creating the data
# y = df["y"]
# x = df.iloc[:, 2:]

# # Building the clustering model and calculating the values of Inertia
# inertias = []
# mapping = {}
# K = range(1, 15)

# for k in K:
#     # Building and fitting the model
#     kmeanModel = KMeans(n_clusters=k).fit(x)
#     kmeanModel.fit(x)
#     inertias.append(kmeanModel.inertia_)
#     mapping[k] = kmeanModel.inertia_

# # Tabulating and Visualizing the Results
# for key, val in mapping.items():
#     print(f'{key} : {val}')

# plt.plot(K, inertias, 'bx-')
# plt.grid(True)
# plt.xlabel('Values of K')
# plt.ylabel('Inertia')
# plt.title('The Elbow Method using Inertia')
# plt.show()

# # Calculate the optimal number of clusters
# y_es =[]
# for i in K:
#     y1 = inertias[len(inertias)-1]
#     y0 = inertias[0]
#     var_y = y1 - y0

#     x1 = K[-1]
#     x0 = K[0]
#     var_x = x1-x0

#     pte = (var_y/var_x)
#     y = pte*(i-2) + y0
#     y_es.append(y)

# distances = []

# for i in range(len(inertias)-1):
#     dist = distance.euclidean(inertias[i],y_es[i])
#     distances.append(dist)

# distances.append(0)

# for i,value in enumerate(distances):
#     if value == max(distances):
#         optimal = i+2

# print('The optimal number of clusters are: ', optimal)

# ##############################################################################
# ################################ CLUSTERING ##################################
# ##############################################################################

# # Create a value for k
# k = optimal

# # Fit and plot the data for each k value
# kmeans = KMeans(n_clusters=k, \
# 				init='k-means++', random_state=42)
# y_kmeans = kmeans.fit_predict(x)

# # Add the ndarray as a new column in the DataFrame
# df['Cluster'] = y_kmeans

# # # MOVE THE LAST COLUMN TO THE 3TH POSITION # #

# # Extract last column
# last_column = df.pop(df.columns[-1])

# # Insert last column in the 4th position
# df.insert(2, last_column.name, last_column)

# print('Saving the file...')
# df.to_csv(PATH+name+"-data_clustered.csv", index=False, sep=';')

# print('The process has been finished :)')

###############################################################################
################################ PATH AND FILE ################################
###############################################################################

df = pd.read_csv(PATH+name+"-data_clustered.csv", sep=";")

###############################################################################
################################# SELECT DATA #################################
###############################################################################

smiles_list = list(df['SMILES'])
molwt_list = []

for sm in smiles_list:
    mol = Chem.MolFromSmiles(sm)
    mol_wt = Chem.Descriptors.MolWt(mol)
    molwt_list.append(mol_wt)

# Add Molecular Weight column
df = df.assign(MolecularWeight=molwt_list)

# Extract last column
last_column = df.pop(df.columns[-1])

# Insert last column in the 2 position
df.insert(3, last_column.name, last_column)

# Sort "MolecularWeight" from smallest to largest
df = df.sort_values(by='MolecularWeight')

# -------------- 0 CLUSTER ------------------ #

# Take only the Cluester 0 values
df_filtered_0 = df[df['Cluster'] == 0]

# Define number of split parts
parts = np.array_split(df_filtered_0, 11)

# Save the results in a new DataFrame
result_df_0 = pd.DataFrame()

# Remove values
for part in parts:

    num_values_to_remove = len(part) - 5

    index_to_remove = np.random.choice(part.index, size=num_values_to_remove, replace=False)

    part = part.drop(index_to_remove)

    result_df_0 = pd.concat([result_df_0, part])

# # -------------- 1 CLUSTER ---------------- #

# # Take only the Cluester 1 values
# df_filtered_1 = df[df['Cluster'] == 1]

# # Define number of split parts
# parts = np.array_split(df_filtered_1, 10)

# # Save the results in a new DataFrame
# result_df_1 = pd.DataFrame()

# # Remove values
# for part in parts:

#     num_values_to_remove = len(part) - 1

#     index_to_remove = np.random.choice(part.index, size=num_values_to_remove, replace=False)

#     part = part.drop(index_to_remove)

#     result_df_1 = pd.concat([result_df_1, part])

# -------------- 2 CLUSTER ---------------- #

# Take only the Cluester 2 values
df_filtered_2 = df[df['Cluster'] == 2]

# Define number of split parts
parts = np.array_split(df_filtered_2, 13)

# Save the results in a new DataFrame
result_df_2 = pd.DataFrame()

# Remove values
for part in parts:

    num_values_to_remove = len(part) - 16

    index_to_remove = np.random.choice(part.index, size=num_values_to_remove, replace=False)

    part = part.drop(index_to_remove)

    result_df_2 = pd.concat([result_df_2, part])

# # -------------- 3 CLUSTER ---------------- #

# # Take only the Cluester 3 values
# df_filtered_3 = df[df['Cluster'] == 3]

# # Define number of split parts
# parts = np.array_split(df_filtered_3, 3)

# # Save the results in a new DataFrame
# result_df_3 = pd.DataFrame()

# # Remove values
# for part in parts:

#     num_values_to_remove = len(part) - 33

#     index_to_remove = np.random.choice(part.index, size=num_values_to_remove, replace=False)

#     part = part.drop(index_to_remove)

#     result_df_3 = pd.concat([result_df_3, part])

# -------------- 4 CLUSTER ---------------- #

# Take only the Cluester 4 values
df_filtered_4 = df[df['Cluster'] == 4]

# Define number of split parts
parts = np.array_split(df_filtered_4, 5)

# Save the results in a new DataFrame
result_df_4 = pd.DataFrame()

# Remove values
for part in parts:

    num_values_to_remove = len(part) - 2

    index_to_remove = np.random.choice(part.index, size=num_values_to_remove, replace=False)

    part = part.drop(index_to_remove)

    result_df_4 = pd.concat([result_df_4, part])

# Filter values
df_filtered_rest = df[df['Cluster'].isin([1,3])]

# Concat dataframes
df_final = pd.concat([result_df_0, df_filtered_rest, result_df_2, result_df_4])

df_final.to_csv(PATH+name+'-data_clustered_filtered.csv', index=False, sep=";")

# ###############################################################################
# ################################ PATH AND FILE ################################
# ###############################################################################

# df_original = pd.read_csv(PATH+name+"-initial_reduction.csv", sep=";")
# df_clustered = pd.read_csv(PATH+name+"-data_clustered_filtered.csv", sep=";")

# ###############################################################################
# ################################# FUSION DATA #################################
# ###############################################################################

# # Take only the y 1 values
# df_original = df_original[df_original['y'] == 1]

# # #Take more than one value
# # df_original = df_original[df_original['y'].isin([1, 2])]

# df_clustered = df_clustered.drop('Cluster', axis=1)
# df_clustered = df_clustered.drop('MolecularWeight', axis=1)

# # Concat dataframes
# df_final = pd.concat([df_clustered, df_original])

# df_final.to_csv(PATH+name+'_final.csv', index=False, sep=";")

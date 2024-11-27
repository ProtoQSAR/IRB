# Load the kmeans model and save the labels in a variable called labels
#%%
import pickle
import numpy as np
import pandas as pd

#%%
# Load the model
kmeans_40k = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/kmeans_model_40k.sav', 'rb'))

labels_40k = kmeans_40k.labels_
counts_labels_40k = pd.Series(labels_40k).value_counts()

len(aa[aa>=3])
morethan3 = aa[aa>=3]
sum(morethan3)

print("for 40k, we have ", len(counts_labels_40k[counts_labels_40k<3]), " clusters with less than 3 compounds")
print("for 40k, we have ", sum(counts_labels_40k[counts_labels_40k<3]), " compounds in clusters with less than 3 compounds")
print("for 40k, we have ", len(counts_labels_40k[counts_labels_40k>=3]), " clusters with 3 or more compounds")
print("for 40k, we have ", sum(counts_labels_40k[counts_labels_40k>=3]), " compounds in clusters with 3 or more compounds")

# %%same for 45k 
kmeans_45k = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/kmeans_model_45k.sav', 'rb'))
labels_45k = kmeans_45k.labels_
counts_labels_45k = pd.Series(labels_45k).value_counts()

print("for 45k, we have ", len(counts_labels_45k[counts_labels_45k<3]), " clusters with less than 3 compounds")
print("for 45k, we have ", sum(counts_labels_45k[counts_labels_45k<3]), " compounds in clusters with less than 3 compounds")
print("for 45k, we have ", len(counts_labels_45k[counts_labels_45k>=3]), " clusters with 3 or more compounds")
print("for 45k, we have ", sum(counts_labels_45k[counts_labels_45k>=3]), " compounds in clusters with 3 or more compounds")

# %%same for 50k
kmeans_50k = pickle.load(open('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/clusters_kmeans/kmeans_model_50k.sav', 'rb'))
labels_50k = kmeans_50k.labels_
counts_labels_50k = pd.Series(labels_50k).value_counts()

print("for 50k, we have ", len(counts_labels_50k[counts_labels_50k<3]), " clusters with less than 3 compounds")
print("for 50k, we have ", sum(counts_labels_50k[counts_labels_50k<3]), " compounds in clusters with less than 3 compounds")
print("for 50k, we have ", len(counts_labels_50k[counts_labels_50k>=3]), " clusters with 3 or more compounds")
print("for 50k, we have ", sum(counts_labels_50k[counts_labels_50k>=3]), " compounds in clusters with 3 or more compounds")

# %%

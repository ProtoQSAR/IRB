# -*- coding: utf-8 -*-
"""
Created on September 2024

@author: Carmen Ortiz - modified proto
usage: python filter_data_chEMBL.py inputpath outputpath_intermediatefile outputpath_finalfile

guide hERG: python filter_data_chEMBL_IRB.py ../datasets_for_modelling/hERG_chEMBL_non_curated.csv ../datasets_for_modelling/intermediate_files/hERGinh_chEMBL_non_curated_indexed.csv ./pre_preprocessed_ONGOING/input_data/TK_ChemBL_hERGinh_firstprocessing.csv

guide Nav1.5: python filter_data_chEMBL_IRB.py ../datasets_for_modelling/Nav_chEMBL_non_curated.csv ../datasets_for_modelling/intermediate_files/Nav15inh_chEMBL_non_curated_indexed.csv ./pre_preprocessed_ONGOING/input_data/TK_ChemBL_Nav15inh_firstprocessing.csv

guide Cav1.5: python filter_data_chEMBL_IRB.py ../datasets_for_modelling/Cav_chEMBL_non-curated.csv ../datasets_for_modelling/intermediate_files/Cav12inh_chEMBL_non_curated_indexed.csv ./pre_preprocessed_ONGOING/input_data/TK_ChemBL_Cav12inh_firstprocessing.csv

ensure that the path exists, as it has been customized for IRB project


"""



import sys
import pandas as pd
#import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def transform_value(row):
    value = row["Standard Value"]
    new_value = value/1000
    row["Standard Value"] = new_value

    return row


def get_class_values(row):
    value = row["Standard Value"]

    if value <= 10:
        new_value = 1
    else:
        new_value = 0

    return new_value

csv_name = str(sys.argv[1])
file_name_no_ext = csv_name.split(".")[0]
# model_name = str(sys.argv[2]) #sligthly modified by eva

outputpath_intermediatefile = str(sys.argv[2])
outputpath_finalfile = str(sys.argv[3])

df = pd.read_csv(csv_name, sep=";")
print(df)

df.index.name = 'ID' #ESC

df.to_csv("%s" % outputpath_intermediatefile, sep = ";", index = True) #ESC index= False --> index = True

#We are going to remove molecules with empty SMILES value.
mask_smiles = pd.notnull(df["Smiles"])
df = df[mask_smiles]

print(df)

#Now we are going to remove empty values from IC50.
mask_ic50 = pd.notnull(df["Standard Value"])
df = df[mask_ic50]

print(df)
#print(df["Standard Units"].value_counts())

df = df.loc[df["Standard Relation"] == "'='"]
print(df)
df = df.loc[df["Standard Units"] == "nM"]
print(df)

print(df["Standard Value"])
#Transform IC50 values from nM to microM
df = df.apply(transform_value, axis = 1)

print(df["Standard Value"])

#Create a new column with classification values for inhibitio --> IC50 <= 10microM = 1 (blocker); IC50 > 10 microM = 0 (non-blocker).
df["y"] = df.apply(get_class_values, axis = 1)

print(df["y"])

df.rename(columns={"Smiles": "SMILES"}, inplace = True)

print(df["SMILES"])

# df[["SMILES", "y"]].to_csv("%s_chEMBL.csv" % model_name, sep = ";", index = False) # modified by proto to rename custom

df[["SMILES", "y"]].to_csv("%s" % outputpath_finalfile, sep = ";", index = True) #ESC index= False --> index = True


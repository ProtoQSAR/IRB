# -*- coding: utf-8 -*-
"""
Created on September 2024

@author: Carmen Ortiz - modified proto
usage: python filter_data_CtoxPred_IRB.py inputpath  outputpath_finalfile

guide Cav12inh: python filter_data_CtoxPred_IRB.py ../datasets_for_modelling/data_cav_dev.csv ./pre_preprocessed_ONGOING/input_data/TK_CtoxPred_Cav12inh_firstprocessing.csv

guide Nav15inh: python filter_data_CtoxPred_IRB.py ../datasets_for_modelling/data_nav_dev.csv ./pre_preprocessed_ONGOING/input_data/TK_CtoxPred_Nav15inh_firstprocessing.csv

ensure that the path exists, as it has been customized for IRB project


"""

import sys
import pandas as pd
import numpy as np

#Function to transform IC50 into a classification dataset. 
def get_class_values(row):
    value = row["pIC50"]

    if value >= 5:
        new_value = 1
    else:
        new_value = 0

    return new_value

#Remove empty values
csv_name = str(sys.argv[1])
file_name_no_ext = csv_name.split(".")[0]
outputpath_finalfile = str(sys.argv[2])

df = pd.read_csv(csv_name, sep=",")
print(df)

#We are going to remove molecules with empty SMILES value. 
mask_smiles = pd.notnull(df["SMILES"])
df = df[mask_smiles]

print(df)

#Now we are going to remove empty values from IC50.
mask_pic50 = pd.notnull(df["pIC50"])
df = df[mask_pic50]

print(df)

df["y"] = df.apply(get_class_values, axis = 1)

print(df["y"].value_counts())

df[["InChl Key", "ChEMBL ID", "SMILES", "y"]].to_csv("%s" % outputpath_finalfile, sep = ";", index = False)

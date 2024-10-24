#%%
"""
Created on October 2024

@author: Carmen. Adapted: proto

"""


############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i) 
##############################################################################

################################ WORK IN SCRIPT PATH ##########################
import os
from pathlib import Path
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')
##############################################################################

################################ WORK IN SCRIPT PATH ##########################
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
#%%
def mw_calculation(row):
    smiles = row["SMILES"]
    m = Chem.MolFromSmiles(smiles)
    mw = Descriptors.ExactMolWt(m)
    row["MW"] = mw
    return row


#%%
file_name = "compare_ONGOING/results/TOX_hERGinh.csv"
df = pd.read_csv(file_name, sep = ";")
print(df)



# %%
df = df.apply(mw_calculation, axis=1)
# %%
print(df["MW"])
# %%
df["MW"].hist(bins=50)
# %%

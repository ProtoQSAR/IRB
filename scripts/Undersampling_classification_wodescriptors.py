#%%
"""
Created on October 2024

@author: Carmen + Ágata --> proto

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
import seaborn as sns


#%%
def mw_calculation(row):
    smiles = row["SMILES"]
    m = Chem.MolFromSmiles(smiles)
    mw = Descriptors.ExactMolWt(m)
    row["MW"] = mw
    return row


#%%
# name = "TOX_Nav15inh"

name = "TOX_hERGinh"

file_name = f"compare_ONGOING/results/{name}.csv"
df = pd.read_csv(file_name, sep = ";")
print(df)



# %%
df = df.apply(mw_calculation, axis=1)
# %%
print(df["MW"])
# %%
df["MW"].hist(bins=50)
# %%


###############################################################################
############################# OUTLIERS ELIMINATION ############################
###############################################################################

unique_values = df['y'].unique()

if set(unique_values) == {0, 1}:
    colors = ['tab:blue', 'tab:red']
elif set(unique_values) == {0, 1, 2}:
    colors = ['tab:blue', 'tab:red', 'tab:green']
else:
    print("Column 'y' has no valid values")

# Calculate Molecular Weight
df['MolecularWeight'] = df['SMILES'].apply(lambda x: Chem.Descriptors.MolWt(Chem.MolFromSmiles(x)))

# Calculate quartiles
q3 = df['MolecularWeight'].quantile(0.75)
q1 = df['MolecularWeight'].quantile(0.25)
iqr = q3 - q1

# Calculate upper and lower limits using interquartile range
upper_limit = q3 + 3 * iqr
lower_limit = q1 - 3 * iqr

# Filter outliers
filtered_outliers = df[(df['MolecularWeight'] > upper_limit) | (df['MolecularWeight'] < lower_limit)]
filtered_outliers.to_csv(f"{name}_outliers.csv", sep=";", index=False)

# Plot box-whisker plot with outliers
sns.set(rc={'figure.figsize':(20,15)})
ax = sns.boxplot(data=df, x="y", y="MolecularWeight", palette=colors)

# Highlight outliers
ax.scatter(filtered_outliers['y'], filtered_outliers['MolecularWeight'], color='red', marker='o', s=100, label='Outliers')

plt.title(f"{name}_highlight_outliers")
plt.legend()
plt.show()

# Remove outliers
remove_outliers = df[(df['MolecularWeight'] >= lower_limit) & (df['MolecularWeight'] <= upper_limit)]
remove_outliers.to_csv(f"{name}_removed_outliers.csv", sep=";", index=False)
df = remove_outliers

# Plot box-whisker plot without outliers
sns.set(rc={'figure.figsize':(20,15)})
sns.boxplot(data=remove_outliers, x="y", y="MolecularWeight", palette=colors).set_title(f"{name}_removed_outliers")

plt.show()

# %%

###############################################################################
################################ RATIO POS/NEG ################################
###############################################################################

print('\t\t negative', df['y'].value_counts()[0])
print('\t\t positive', df['y'].value_counts()[1])





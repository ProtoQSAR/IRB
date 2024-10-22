"""
Description: A Script to perform BitBirch clustering on a dataset of molecules using MACCS fingerprints.
It reads a CSV file with SMILES strings, calculates MACCS fingerprints, and saves the fingerprints in a pickle file.

Author: Laureano E. Carpio (Moldrug)

Date: 2024-10-21

"""

#%% imports
import pandas as pd
import numpy as np

from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import AllChem

import sys
sys.path.insert(1, '/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/BitBirch/')
import bitbirch

# %%
bb = bitbirch.BitBirch(threshold=0.95)

fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3)
# %%
def smi2numpyarr(smi):
    mol = Chem.MolFromSmiles(smi)
    fp = fpgen.GetFingerprintAsNumPy(mol)
    return fp

df = pd.read_csv("/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/clustering/BitBirch/IRB_library_merged_all-prepared.csv", sep = ";")

fpnp = np.array([smi2numpyarr(smi) for smi in df.SMILES.to_list()])

fpnp.shape

# %%
%time res = bb.fit(fpnp)
# %%
num_clust_g3 = sum(1 for c in res.get_cluster_mol_ids() if len(c) > 3)

print("# bb clusters with >3 compounds: ", num_clust_g3)
# %%
centroids = res.get_centroids()
# %%
mol_ids = res.get_cluster_mol_ids()
# %%
mol_ids
# %%
print(mol_ids)
# %%

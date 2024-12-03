# -*- coding: utf-8 -*-
"""
Created on Fri May 12 06:37:23 2023

@author: dataco
"""

#%%
####input formatting

import pandas as pd

PPARg = pd.read_excel(r'C:\Users\dataco\protoqsar\to_do\Para David\validation\PPARg.xlsx', sheet_name =['ChEMBL', 'toxCSM'])
AhR = pd.read_excel(r'C:\Users\dataco\protoqsar\to_do\Para David\validation\AhR.xlsx', sheet_name =['ChEMBL', 'toxCSM'])
NRF2 = pd.read_excel(r'C:\Users\dataco\protoqsar\to_do\Para David\validation\NRF2.xlsx')
# Mitotox
PXR = pd.read_excel(r'C:\Users\dataco\protoqsar\to_do\Para David\validation\PXR.xlsx')
PATH = "C:/Users/dataco/protoqsar/to_do/Para David/validation/"
INPUT = ['PPARg.xlsx','AhR.xlsx', 'NRF2.xlsx', 'PXR.xlsx' ]
DATAS = PPARg,AhR,NRF2,PXR
NAMES = []
OUTPUT=[]
OUTPUTFOLDER1 = 'noheader/'
OUTPUTFOLDER2 = 'header/'
OUTPUTFOLDER3 = 'fasta/'
for nam in INPUT:
    NAMES.append(nam.split('.')[0])


#%%
#fasta

#AhR
#split = 500

ds1=PPARg.get('ChEMBL')
ds2=PPARg.get('toxCSM')

ds1 = ds1['SMILES']
ds2 = ds2['SMILES']


arc = open("smiles_PPARg_ChEMBL.fasta", "x")
for i,smi in enumerate(ds1):
    arc.write('>SMI_{}\n'.format(i))
    arc.write(smi+'\n')
   
arc.close()


#%%


arc = open("smiles_PPARg_toxCSM.fasta", "x")
for i,smi in enumerate(ds2):
    arc.write('>SMI_{}\n'.format(i))
    arc.write(smi+'\n')
   
arc.close()




#%%
data = NRF2.SMILES

arc = open("smiles_NRF2.fasta", "x")
for i,smi in enumerate(data):
    arc.write('>SMI_{}\n'.format(i))
    arc.write(smi+'\n')
   
arc.close()


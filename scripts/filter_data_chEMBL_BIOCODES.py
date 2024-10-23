# -*- coding: utf-8 -*-
"""
Created on September 2024

@author: Carmen Ortiz - modified proto
usage: python filter_data_chEMBL.py inputpath outputpath_intermediatefile outputpath_finalfile



guide Cav1.5: python filter_data_chEMBL_BIOCODES.py ../datasets_for_modelling/Abaumannii_MIC_CHEMBL_filter-nM.csv ../datasets_for_modelling/intermediate_files/AbaumanniiMIC_chEMBL_non_curated_indexed.csv ./pre_preprocessed_ONGOING/input_data/BIO_ChemBL_AbaumanniiMIC_firstprocessing.csv

ensure that the path exists, as it has been customized for BIOCODES project


"""



import sys
import pandas as pd
#import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors



csv_name = str(sys.argv[1])
file_name_no_ext = csv_name.split(".")[0]
# model_name = str(sys.argv[2]) #sligthly modified by eva

outputpath_intermediatefile = str(sys.argv[2])
outputpath_finalfile = str(sys.argv[3])

df = pd.read_csv(csv_name, sep=";")


df.index.name = 'ID' #ESC

df.rename(columns = {'Standard Value' : 'y', 'Smiles': 'SMILES'}, inplace = True)



df.to_csv("%s" % outputpath_intermediatefile, sep = ";", index = True) #ESC index= False --> index = True

df[["SMILES", "Molecule ChEMBL ID", "y"]].to_csv("%s" % outputpath_finalfile, sep = ";", index = True) #ESC index= False --> index = True


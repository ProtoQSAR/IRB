#%% A script that creates the folder for_screening and then load a csv file, and split this file in dataframes of 10000 rows, and save the different df as csv separated by ;

import os
import pandas as pd
import numpy as np


def split_df(df, splits,path):
    for i,df in enumerate(np.array_split(df, splits)):
        chunk_name = f"IRB__screening__chunk-{str(i)}"
        print(chunk_name, df.shape)
        outfile = chunk_name+'.csv'
        del df["y"]
        df.to_csv(path+outfile, sep = ';', index = False)
        print(f"File {outfile} saved")

# Create folder for_screening
# create_folder('for_screening')
path = '/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/libreria/finals/for_screening/'

os.makedirs(path)
#%%

# Load csv file

df = pd.read_csv('/home/lauri/Desktop/ProtoQSAR/PROJECTS/contratas/IRB/libreria/finals/IRB_library_merged_all.csv', sep = ';')

# Split the df in 10000 rows chunks
split_df(df, 14, path)

# %%
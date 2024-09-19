# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:13:33 2024

@author: proto
"""

import os

os.chdir(r'C:\Users\Enrique\Documents\GitHub\IRB\datasets_for_modelling')
import pandas as pd


properties = ['cyp1a2_inhibitor','cyp2C9_inhibitor', 'cyp2C19_inhibitor', 'cyp2d6_inhibitor', 'cyp3a4_inhibitor']

for properti in properties:
    df_train = pd.read_csv(f'{properti}_train.csv', sep = ',')
    df_test = pd.read_csv(f'{properti}_test.csv', sep = ',')
    df_val = pd.read_csv(f'{properti}_val.csv', sep = ',')
    
    df_all = pd.concat([df_train,df_test,df_val], axis = 0)
    
    print(properti)
    print('\t\t negative', df_all['label'].value_counts()[0])
    print('\t\t positive', df_all['label'].value_counts()[1])
    
    df_all.index.name = 'ID'
    
    df_all.to_csv(f'DeppPK_{properti}_all.csv', sep = ';', index= True)
import os
import pandas as pd
import sys

file_name = sys.argv[1]

df_train = pd.read_csv(f'{file_name}_train.csv', sep = ',')
df_test = pd.read_csv(f'{file_name}_test.csv', sep = ',')
df_val = pd.read_csv(f'{file_name}_val.csv', sep = ',')

df_all = pd.concat([df_train,df_test,df_val], axis = 0)

print('\t\t negative', df_all['label'].value_counts()[0])
print('\t\t positive', df_all['label'].value_counts()[1])
    
df_all.index.name = 'ID'

df_all.to_csv(f'DeepPK_{file_name}_all.csv', sep = ';', index= True)

import pandas as pd
import os

PATH = "C:/Users/Enrique/Documents/GitHub/IRB/Models/CYP3A4 Inhibitor/calc/"

df = pd.DataFrame()

for root, dirs, files in os.walk(PATH):
    for file in files:
        if "paralel_calculated_with_y" in file:
            file_path = os.path.join(root, file) 
            info = pd.read_csv(file_path, sep=';')
            df = pd.concat([df, info], ignore_index=True) 


output_path = os.path.join(PATH, 'CYP2D6_Inhibitor-paralel_calculated_with_y.csv')
df.to_csv(output_path, sep=';', index=False)

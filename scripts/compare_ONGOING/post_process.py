import pandas as pd
import ntpath

#Read merged file. 
path_file = "/home/carmenortiz/Escritorio/contratas/IRB_repo/scripts/compare_ONGOING/results/all_merged_IRB_hERGinh.csv"
df = pd.read_csv(path_file, sep=";", encoding="utf-8")

#Get the model name.
file_name = ntpath.basename(path_file)
model_name = file_name.split(".")[0].replace("all_merged_", "")
print(model_name)
print(df.shape)

#Remove outliers.
df=df[df["valid/outlier"]=="valid"]

#Get SMILES and endpoint columns.
df=df[["SMILES","y_consensus"]]
print(df.columns)

#Rename columns for WOTAN.
df.rename(columns={'y_consensus': 'y'},inplace=True)
print(df.columns)

#Save file.
df.to_csv("results_processed/%s.csv" % model_name, sep=";", encoding="utf-8", index=False)


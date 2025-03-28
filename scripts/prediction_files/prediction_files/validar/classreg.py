import pandas as pd


############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i)
##############################################################################
#################################### IMPORTS #################################
import os
import re   #"""importar regular expresions"""
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score, confusion_matrix

##############################################################################
################################ INITIAL VARIABLES ###########################

df=pd.read_csv("D:\\GitHub\\IRB\\scripts\\prediction_files\\prediction_files\\validar\\HHcarcino_outputVEGA_allcolumns.csv",sep=";",encoding="utf8")
print(df)

y_col= "Experimental_y"


regression_columns = [col for col in df.columns if col.split("_",1)[0]=="reg"]
classification_columns = [col for col in df.columns if col.split("_",1)[0]=="class"]

df_pair_reg = [df[[y_col, col]] for col in regression_columns]
df_pair_class = [df[[y_col, col]] for col in classification_columns]

print(df_pair_reg[0])
print(df_pair_class[1])


for reg_col in regression_columns:
    experimental_values=df[y_col]
    predicted_values = df[reg_col]

    mse = mean_squared_error(experimental_values, predicted_values)
    r2 = r2_score(experimental_values, predicted_values)
    plt.figure(figsize=(8, 6))
    plt.scatter(experimental_values, predicted_values, alpha=0.6, label="Datos")
    plt.plot(
        experimental_values,
        experimental_values,
        color="red",
        linestyle="--",
        label="y = x",
    )
    plt.text(
        0.05,
        0.80,
        f"$R^2$ = {r2:.2f}",
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.5),
    )
    plt.xlabel("Experimental Values")
    plt.ylabel("Predicted Values")
    plt.title(f"Experimental vs. {col}")
    plt.legend()

    # Normalizar el nombre del archivo reemplazando caracteres inválidos
    safe_col_name = col.replace(":", "_").replace("/", "_")
    plt.savefig(f"correlation_{safe_col_name}.png")
    plt.show()

    # Métricas de ajuste del modelo
    print(f"Column: {col}")
    print(f"Mean Squared Error: {mse:.2f}")
    print(f"R^2 Score: {r2:.2f}")
    print("-" * 40)

for class_col in classification_columns:
    experimental_values=df[y_col]
    predicted_values = df[class_col]


    tn_test, fp_test, fn_test, tp_test = confusion_matrix(
    experimental_values,
    predicted_values,).ravel()
   
    
    y_exp_tag = "Experimental"
    y_pred_tag = "Prediction"
    positive_tag = "Toxic"
    negative_tag = "No_toxic"
    # Imprimir matriz de confusión
print(f"\n\n{y_exp_tag} vs {y_pred_tag}")
print("\n\t\t\t        Predicted")
print(f"\t\t\t\t{negative_tag}\t{positive_tag}")
print(f"Experimental\t{negative_tag}\t" + str(tn_test) + "\t\t" + str(fp_test))
print(f"\t\t{positive_tag}\t" + str(fn_test) + "\t\t" + str(tp_test))

# Calcular métricas de clasificación
accuracy = (tp_test + tn_test) / (tp_test + tn_test + fp_test + fn_test)
sensitivity = tp_test / (tp_test + fn_test)
specificity = tn_test / (tn_test + fp_test)
precision = tp_test / (tp_test + fp_test)
f1 = (2 * sensitivity * precision) / (sensitivity + precision)

accuracy, sensitivity, specificity, precision, f1 = (
    round(accuracy, 2),
    round(sensitivity, 2),
    round(specificity, 2),
    round(precision, 2),
    round(f1, 2),
)

# Imprimir métricas
print("\nAccuracy:", accuracy)
print("Sensitivity:", sensitivity)
print("Specificity:", specificity)
print("Precision:", precision)
print("F1 Score:", f1)

# Métricas adicionales
sensit_TPR_test = tp_test / (tp_test + fn_test)
spec_TNR_test = tn_test / (tn_test + fp_test)
NPV_test = tn_test / (tn_test + fn_test)
FNR_test = fn_test / (fn_test + tp_test)
FPR_test = fp_test / (fp_test + tn_test)
FDR_test = fp_test / (fp_test + tp_test)
FOR_test = fn_test / (fn_test + tn_test)
F_score_test = 2 * tp_test / (2 * tp_test + fp_test + fn_test)
CSI_test = tp_test / (tp_test + fn_test + fp_test)

# Imprimir métricas adicionales
print("\n\t        | Test |")
print("       acc\t| %3.2f |" % (accuracy))
print("  prec_PPV\t| %3.2f |" % (precision))
print("    recall\t| %3.2f |" % (sensit_TPR_test))
print("        f1\t| %3.2f |" % (f1))
print("sensit_TPR\t| %3.2f |" % (sensit_TPR_test))
print("  spec_TNR\t| %3.2f |" % (spec_TNR_test))
print("       NPV\t| %3.2f |" % (NPV_test))
print("       FNR\t| %3.2f |" % (FNR_test))
print("       FPR\t| %3.2f |" % (FPR_test))
print("       FDR\t| %3.2f |" % (FDR_test))
print("       FOR\t| %3.2f |" % (FOR_test))
print("   F_score\t| %3.2f |" % (F_score_test))
print("       CSI\t| %3.2f |" % (CSI_test))

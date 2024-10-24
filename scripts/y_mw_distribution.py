#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
#%%
def mw_calculation(row):
    smiles = row["SMILES"]
    m = Chem.MolFromSmiles(smiles)
    mw = Descriptors.ExactMolWt(m)
    row["MW"] = mw
    return row


#%%
file_name = "TOX_MRDD-paralel_calculated_with_y.csv"
df = pd.read_csv(file_name, sep = ";")
print(df)
# %%
df["y"].hist(bins=50)
# %%
y_values = df['y']  # Asumiendo que la columna se llama 'y'
y_values_positive = y_values + 1  # Para evitar valores negativos o cero
# y_transformed_boxcox, lambda_boxcox = boxcox(y_values_positive)
# Aplicar Transformación Logarítmica
y_transformed_log = np.log10(y_values)
y_transformed_neglog=np.log10(1/y_values)
# Aplicar Transformación de Raíz Cuadrada
y_transformed_sqrt = np.sqrt(y_values_positive)
# Aplicar Transformación de Raíz Cuadrada Inversa
y_transformed_sqrt_inv = np.sqrt(np.sqrt(y_values_positive))
# Aplicar Transformación Recíproca
y_transformed_reciprocal = 1 / (y_values_positive)
# Aplicar Transformación de Yeo-Johnson (funciona con valores negativos y positivos)
# y_transformed_yeojohnson, lambda_yeojohnson = yeojohnson(y_values)
# Visualizar todas las transformaciones
plt.figure(figsize=(18, 12))
plt.subplot(3, 2, 1)
plt.hist(y_values, bins=100, color='blue', alpha=0.7)
plt.title('Distribución Original de y')
plt.xlabel('Valores de y')
plt.ylabel('Frecuencia')
plt.subplot(3, 2, 2)
plt.hist(y_transformed_neglog, bins=100, color='green', alpha=0.7)
plt.title('Transformación de y log 1/y')
plt.xlabel('')
plt.subplot(3, 2, 3)
plt.hist(y_transformed_log, bins=100, color='purple', alpha=0.7)
plt.title('Transformación de y (Logarítmica)')
plt.xlabel('Valores Transformados de y (Log)')
plt.subplot(3, 2, 4)
plt.hist(y_transformed_sqrt, bins=100, color='orange', alpha=0.7)
plt.title('Transformación de y (Raíz Cuadrada)')
plt.xlabel('Valores Transformados de y (Raíz Cuadrada)')
plt.subplot(3, 2, 5)
plt.hist(y_transformed_sqrt_inv, bins=100, color='red', alpha=0.7)
plt.title('Transformación de y (Raíz Cuadrada Inversa)')
plt.xlabel('Valores Transformados de y (Raíz Cuadrada Inversa)')
plt.subplot(3, 2, 6)
plt.hist(y_transformed_reciprocal, bins=100, color='brown', alpha=0.7)
plt.title('Transformación de y (Recíproca)')
plt.xlabel('Valores Transformados de y (Recíproca)')
plt.tight_layout()
plt.show()
# %%
df = df.apply(mw_calculation, axis=1)
# %%
print(df["MW"])
# %%
df["MW"].hist(bins=50)
# %%

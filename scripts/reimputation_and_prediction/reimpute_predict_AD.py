# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:37:05 2022

@author: Rita. Adapted: proto

@author: proto
"""


############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i) 
##############################################################################


################################ WORK IN SCRIPT PATH ##########################
import os
from pathlib import Path

##############################################################################

###############################################################################
##################################IMPORTS######################################
###############################################################################

import pandas as pd
import numpy as np

import pickle
from sklearn import set_config
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler, MinMaxScaler

########################## IMPORTS FOR MODELS##################################

from sklearn.tree import DecisionTreeClassifier
from lightgbm import LGBMRegressor

from sklearn.svm import SVR, LinearSVR, NuSVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import LinearRegression, Lasso, SGDRegressor, Ridge, PassiveAggressiveRegressor, LassoLars
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    ExtraTreesRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor
)
from sklearn.neural_network import MLPRegressor
# import xgboost as xgb
import lightgbm as lgb

# classification algorithms
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
# from sklearn.naive_bayes import GaussianNB, MultinomialNB # Not working with negative values
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression, Lasso, SGDClassifier
from sklearn.neural_network import MLPClassifier
import lightgbm as lgb
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import (
    RandomForestClassifier,
    ExtraTreesClassifier,
    AdaBoostClassifier,
    GradientBoostingClassifier
)

########################  IMPORTS FOR METRICS##################################

# from sklearn.metrics import r2_score
# from sklearn.metrics import mean_squared_error
# from sklearn.metrics import mean_absolute_error
# from sklearn.metrics import explained_variance_score

# import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
import json
from rdkit import Chem
from utils.get_metrics import getMetrics
########################  IMPORTS FOR  PARALLELIZATION ########################


# from AD import calculate_ADs
import multiprocessing
from multiprocessing import Process
from descriptors import wotan
import time
##############################################################################
##############################################################################
##############################################################################

################################ INITIAL VARIABLES ###########################
input_data_folder =  '..' + os.path.sep + '..' + os.path.sep + 'Models'

##############################################################################

#%% FUNCTIONS

class Calculate(Process):
    def __init__ (self, input_df,i,ini,fin,model_dict,return_dict):
        

        Process.__init__ (self)
        self.i = i
        self.input_df = input_df
        self.ini = ini
        self.fin = fin
        self.model_dict = model_dict
        self.return_dict = return_dict



    def run (self):

        self.df = self.input_df.iloc[self.ini:self.fin,:]

        molec_struc = [Chem.MolFromSmiles(smi) for smi in self.df['SMILES']]
        smiles = [smi for smi in self.df['SMILES']]


        self.results = pd.DataFrame({'SMILES': smiles})
        
        # self.resultsresults = pd.DataFrame()

        for desc in self.model_dict:
            kwargs = self.model_dict[desc]
    
            desc_fn_name = desc.split('$')[0]
            desc_name = desc.split('$')[1]
            desc_fn = getattr(wotan, desc_fn_name)
    
            lista = []
            for mol in molec_struc:
    
    
                try:
                    res_per_desc = desc_fn(mol, **kwargs)
                    lista.append(res_per_desc)
    
                except:
                    lista.append('NaN')
    
            self.results[desc_name] = pd.Series(lista)



        self.return_dict[self.i] = self.results





### in case that non-parallel option is desired
def calculate_all(molec_struc, model_dict, wotan_path):



    results = pd.DataFrame()

    for desc in model_dict:

        kwargs = model_dict[desc]

        desc_fn_name = desc.split('$')[0]
        desc_name = desc.split('$')[1]
        desc_fn = getattr(wotan, desc_fn_name)

        lista = []
        for mol in molec_struc:


            try:
                res_per_desc = desc_fn(mol, **kwargs)
                lista.append(res_per_desc)

            except:
                lista.append('NaN')

        results[desc_name] = pd.Series(lista)
        # results[desc_name] = [desc_fn(mol, **kwargs) for mol in molec_struc]

#     return results


def load_json(json_path):
    '''Loads JSON files as a dictionary.
    '''
    with open(json_path, 'r', encoding='utf8') as json_file:
        json_dict = json.load(json_file)

    return json_dict

def nanify(feats_df):
    ### substitute text in descriptors by NaN
    print("antes de nanify", feats_df.shape)

    df_nan = feats_df.replace(to_replace=r"[a-zA-Z]",
                                value=np.nan,
                                regex=True).replace(np.inf,np.nan)
    to_cast=list(df_nan.select_dtypes(include=[object]).columns)
    df_nan[to_cast] = df_nan[to_cast].astype(dtype=np.float64)
    df_nan[df_nan <= -1e38] = np.nan
    df_nan[df_nan >= 1e38] = np.nan

    print("después de nanify", df_nan.shape)

    return df_nan

def Knn_imputation(dataset):

    # This is a Copied version of the Knn_imputation function in wotan
    # The only difference (I think) is that does not remove the first columns, this is done before
    from time import perf_counter

    set_config(working_memory=4)

    # print(dataset.shape)

    features = dataset

    out_columns = list(features.columns)  # SME EDIT

    df_nan = features.replace(
        to_replace=r'[a-zA-Z]', value=np.nan, regex=True).replace(np.inf, np.nan)

    to_cast = list(df_nan.select_dtypes(include=[object]).columns)
    df_nan[to_cast] = df_nan[to_cast].astype(dtype=np.float64)

    df_nan[df_nan <= -1e38] = np.nan
    df_nan[df_nan >= 1e38] = np.nan

    all_nan = df_nan.columns[df_nan.isna().all()].tolist()
    if len(all_nan) > 0:  # SME EDIT
        print('Warning! Some of descriptors are NaN along the full data set:', len(
            all_nan), 'BIC0_core' in all_nan)  # SME EDIT
    for i in all_nan:  # SME EDIT
        out_columns.remove(i)  # SME EDIT
    df_nan = df_nan.dropna(axis=1, how='all')

    #df_nan = df_nan.round(3)

    t1 = perf_counter()
    imputer = KNNImputer(missing_values=np.nan,
                         n_neighbors=3, weights="uniform")
    print('[+] fitting')
    result = imputer.fit(df_nan)
    print('[+] transforming')
    result = imputer.transform(df_nan)

    # print(result.shape,len(list(features.columns)))  #SME EDIT
    imputed_df = pd.DataFrame(result, columns=out_columns)  # SME EDIT
    t2 = perf_counter()
    # print((t2-t1)/60)
    set_config(working_memory=1024)
    return imputed_df, imputer




def calc_descriptors(option, mod_name, input_df, model_dict, model, set_data, output_path):
    
    

    if option == "reffiting":

        ######### in case that non-parallel option is desired #################
        # print('[+] Calculating w/o parallelization')
        # start_total_time = time.time()
        # mols = [Chem.MolFromSmiles(mol) for mol in list(validation_data['SMILES'])]
        
        # descriptors = calculate_all(mols, model_dict, wotan_path)
        # # print(descriptors)
        
        # elapsed_time = time.time() - start_total_time
        # print(f'{round(elapsed_time, 2)} seconds')
        #######################################################################
        
        start_total_time = time.time()
        print('\t\t[+++] Calculating descriptors with parallelization')

        cores = multiprocessing.cpu_count()

        NUM_PROCESOS = 3*int(cores/4)
        # NUM_PROCESOS = 1
        shape = input_df.shape
        
        interv = int(shape[0]/NUM_PROCESOS)

        jobs = []

        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        for i in range(NUM_PROCESOS):

            if i == NUM_PROCESOS-1:
                ini =  i*interv
                fin = shape[0]
            else:
                ini =  i*interv
                fin = (i+1)*interv

            jobs.append(Calculate(input_df,i,ini,fin,model_dict,return_dict))
            jobs[i].start()

        for job in jobs:
            job.join()


        descriptors_df = pd.DataFrame()

        for i in range(NUM_PROCESOS):
            descriptors_df = pd.concat([descriptors_df,return_dict[i]], axis=0, sort=False)

        #print(final_df,data_df)
        # final_df = pd.concat([final_df.reset_index(drop=True),data_df[extra_data]], axis=1, sort=False) 
        
        # print(descriptors_df)
        elapsed_time = time.time() - start_total_time
        print(f'{round(elapsed_time, 2)} seconds')


        descriptors_df.to_csv(output_path + os.path.sep + f'raw_{mod_name}-descriptors-train.txt', sep = '\t', index = False)

        return descriptors_df


def refit ():
    print('pppp')


#%%



if __name__ == '__main__':
    
    parent = Path(__file__).resolve().parent
    os.chdir(parent)
    print(f'Working on: \n {os.getcwd()}')
    
    df_models_to_process = pd.read_excel('../reimputation_and_prediction_data.xlsx')
    
    print('\n[!] Yo are in reimputation mode[!]')
    
    option = 'reffiting'
    '''
    reffiting --> reffits locally a previus model
                     needs training-set, test-set, sav file, json file
    prediction --> predicts a new dataset
                   needs training-set, test-set, sav file, json file
    '''
    df_metrics = pd.DataFrame()
    ### IF mode = = reimpute
    if option == 'reffiting':
     
        df_models_to_process = df_models_to_process[df_models_to_process['reffit'] == 'yes']
        

        
    for mod_name in df_models_to_process['final_model_name']:
        
        print(f'\n[+] Working on {mod_name}')
        
        print('\t[++] Reading input files')
        #read input files
        
        path_model = input_data_folder + os.path.sep + mod_name
        
        train_df = pd.read_csv(path_model + os.path.sep + f'{mod_name}-descriptors-train.txt', sep = '\t')
        test_df = pd.read_csv(path_model + os.path.sep + f'{mod_name}-descriptors-test.txt', sep = '\t')
        model = pickle.load(open(path_model + os.path.sep + f'{mod_name}.sav', "rb"))
        model_dict =  load_json(path_model + os.path.sep + f'{mod_name}.json')
        scaler =  pickle.load(open(path_model + os.path.sep + f'{mod_name}-train_set.sca', "rb"))
       
        print(f'\t\t Model type: {model._estimator_type}')
        
        print('\t[++] Calculating original metrics')
        
        df_metrics_byone = getMetrics(mod_name, model, train_df, test_df, option)
        
        df_metrics = pd.concat([df_metrics,df_metrics_byone], axis = 0)


        print('\t[++] Refitting model')

        results_folder = path_model + os.path.sep + f'imputed-{mod_name}'

        if not os.path.exists(results_folder):
            os.makedirs(results_folder)
        
        train_descriptors = calc_descriptors(option, mod_name, train_df, model_dict, model, 'train', results_folder)
        test_descriptors = calc_descriptors(option, mod_name, test_df, model_dict, model, 'test', results_folder)
            


        print('\t\t[+++] Imputing')
        
        descriptors = train_descriptors.iloc[:,1:]
        
    
        imputed_df, imputer = Knn_imputation(descriptors)
    
        pickle.dump(imputer, open(output_path + os.path.sep + f'{mod_name}-train_set.imp', 'wb'))
        
        
        
        print('\t\t[+++] Scaling')
        
        
    
        
        print('\t\t[+++] Fitting new model')
 
    # #Scaler

    # sc_file = f"{mod_name}-train_set.sca"
    # old_scaler = pickle.load(open(PATH + sc_file,"rb"))
    # if isinstance(old_scaler, StandardScaler):
    #     print(f"cojo el {mod_name}-train_set.sca")
    #     scaler = StandardScaler()
    # elif isinstance(old_scaler, MinMaxScaler):
    #     scaler = MinMaxScaler()
    # else:
    #     print(f"Scaler file {sc_file} not found.")
    #     print("StandardScaler used by default")
    #     scaler = StandardScaler()
    # # scaler = StandardScaler()
    # train2scale = pd.read_csv(PATH+f"{mod_name}-train_set.csv", sep=";")
    # train2scale_feats = train2scale[sel_feats]
    # train_scaled = scaler.fit(train2scale_feats)

    # #impute and scale validation data
    # val_nan = nanify(descriptors)
    # val_feats_imput = imputer.transform(val_nan)
    # val_feats_scaled = scaler.transform(val_feats_imput)

    # if class_reg == 1:
    #     print("classification model")

    #     #predict classes
    #     predictions = model.predict(val_feats_scaled)
    #     pred_prob = model.predict_proba(val_feats_scaled)
    #     pred_prob_max = [max(probs) for probs in pred_prob]

    #     # classification results

    #     val_feats_scaled_df = pd.DataFrame(val_feats_scaled, columns=sel_feats)
    #     results = pd.concat([validation_data.loc[:,["SMILES", "predicted"]], val_feats_scaled_df],axis=1)
    #     results["class_pred"] = predictions
    #     results["class_pred_prob"] = pred_prob_max

    # else:
    #     print("regression model")

    #     #predict pKa
    #     predictions = model.predict(val_feats_scaled)

    #     #regression results
    #     val_feats_scaled_df = pd.DataFrame(val_feats_scaled, columns=sel_feats)
    #     results = pd.concat([validation_data.loc[:,["SMILES", "y"]], val_feats_scaled_df],axis=1)
    #     #results.rename(columns={"y":"pKa_exp"}, inplace=True)
    #     results["y_pred"] = predictions

    # if out_type == "1":
    #     return results
    # elif  out_type == "2":
    #     return scaler, imputer






        
        #set outputs
        
        results_folder = path_model + os.path.sep + f'imputed-{mod_name}'

        if not os.path.exists(results_folder):
            os.makedirs(results_folder)
            print(f"Output directory created: {results_folder}")
        else:
            print(f"Output directory already exists: {results_folder}")
            
            
            
            
            
            
            
    output_metrics_filename =   f'metrics_{time.strftime("%Y%m%d_%H%M%S")}.csv'  
    print(f'\t\t metrics file created: {output_metrics_filename}')
    df_metrics.to_csv(output_metrics_filename, sep = ';')    
            
            
        #obtain metrics original model

        

#%% ADs
# PATH2 = r'C:\Users\proto\Desktop\tothexinxol\IRB\Models\TK_OATP1B1inh_no3D_us\predictions'    
    
  
# train_df = pd.read_csv(f"{PATH2}/{mod_name}-descriptors-train.txt",sep="\t")    
# test_df = pd.read_csv(f"{PATH2}/predict_{mod_name}.csv",sep=";")  
# ads = calculate_ADs(train_df,test_df,mod_name)

# ads.to_csv(f"{PATH2}/predict_{mod_name}_ad.csv",sep=";",index=False)

'''
elapsed_time = 11.042634963989258 30 mols
'''
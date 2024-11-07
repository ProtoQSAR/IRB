# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:37:05 2022

@author: Rita. Adapted: proto


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
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
########################  IMPORTS FOR  PARALLELIZATION ########################


# from AD import calculate_ADs
import multiprocessing
from multiprocessing import Process
from descriptors import wotan
from utils.get_metrics import getMetrics
from utils.AD import calculate_ADs
import time
##############################################################################
##############################################################################
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

    return results

def inputoutput_files(model, path_model):
    
    train_df = pd.read_csv(path_model + os.path.sep + f'{mod_name}-descriptors-train.txt', sep = '\t')
    test_df = pd.read_csv(path_model + os.path.sep + f'{mod_name}-descriptors-test.txt', sep = '\t')
    model = pickle.load(open(path_model + os.path.sep + f'{mod_name}.sav', "rb"))
    model_dict =  load_json(path_model + os.path.sep + f'{mod_name}.json')
    scaler =  pickle.load(open(path_model + os.path.sep + f'{mod_name}-train_set.sca', "rb"))

    return  train_df, test_df, model, model_dict, scaler  


def load_json(json_path):
    '''Loads JSON files as a dictionary.
    '''
    with open(json_path, 'r', encoding='utf8') as json_file:
        json_dict = json.load(json_file)

    return json_dict

def nanify(df):
    ### substitute text in descriptors by NaN


    df_nan = df.replace(to_replace=r"[a-zA-Z]",
                                value=np.nan,
                                regex=True).replace(np.inf,np.nan)
    to_cast=list(df_nan.select_dtypes(include=[object]).columns)
    df_nan[to_cast] = df_nan[to_cast].astype(dtype=np.float64)
    df_nan[df_nan <= -1e38] = np.nan
    df_nan[df_nan >= 1e38] = np.nan



    return df_nan







def calc_descriptors(option, mod_name, input_df, model_dict, set_data, output_path):
    
    



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
    print(f'\t\t[+++] Calculating descriptors with parallelization for {set_data}')

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

    
    elapsed_time = time.time() - start_total_time
    print(f'\t\t\t{round(elapsed_time, 2)} seconds')

    raw_descriptors_file_name = output_path + os.path.sep + f'raw_{mod_name}-descriptors-{set_data}.txt'
    descriptors_df.to_csv(raw_descriptors_file_name, sep = '\t', index = False)
    print(f'\t\t file with raw descriptors created: {raw_descriptors_file_name}')

    return descriptors_df

def impute(df, imputer):
        
    imputed_df = imputer.transform(df)
    
    return imputed_df

def scale(df, scaler):
    
    scaled_df = scaler.transform(df)
    
    return scaled_df
    
def predict(df, model):
    Y_pred = model.predict(df)
    
    return Y_pred

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"\t\tOutput directory created: {path}")
    else:
        print(f"\t\tOutput directory already exists: {path}")    
#%%



if __name__ == '__main__':
    
    parent = Path(__file__).resolve().parent
    os.chdir(parent)
    print(f'Working on: \n {os.getcwd()}')
    
    
    
    
    
    option = 'reffit'
    
    calculate = 'WITH'
    
    '''
    
    all of them need reimputation_and_prediction_data.xlsx file
    reffit --> reffits locally a previus model
                     needs training_df, test_df, sav file, json file, optional: scaler file
                     
                     calculate = 'WITHOUT' --> I have already recalculate new descriptors for train and test
                     calculate = 'WITH' --> I need to calculate them
    calculate --> calculate descriptors of a custom json file 
                   needs new_dataset, json file                     
                     
    predict --> predicts a new dataset with precalculated descriptors
                   needs new_dataset, training_df, test_df, sav file, json file, imp file, scaler file 
                   
                   
    mergeJSONS --> create a consensus json file for different models 
                   needs a set of json files
                   
    '''

    
    
    df_models_to_process = pd.read_excel('reimputation_and_prediction_data.xlsx', sheet_name = 'models')

    df_models_to_process = df_models_to_process[df_models_to_process[option] == 'yes']    
    

    if option != 'reffit' :
        
        
        print(f'\n[!] You are in {option} mode [!] ')
        
        
        if option != 'mergeJSONS':
            
            df_df_to_process = pd.read_excel('reimputation_and_prediction_data.xlsx', sheet_name = 'df_to_predict')
            
            df_df_to_process = df_df_to_process[df_df_to_process[option] == 'yes']
        
        
        
    else :
        print(f'\n[!] You are in {option} mode {calculate} calculating descriptors [!]')
    
        
    ##########################################################################     
        
        
    if option == 'reffit':
        
        df_metrics = pd.DataFrame() 
        
        for mod_name in df_models_to_process['final_model_name']:
            
            print(f'\n[+] Working on {mod_name} model')
            
            print('\t[++] Reading input files and creating putput directory')

            model_data_folder =  '..' + os.path.sep + '..' + os.path.sep + 'Models'
            path_model = model_data_folder + os.path.sep + mod_name
            results_folder = path_model + os.path.sep + f'imputed-{mod_name}'
    
            makedir(results_folder)
                
                
            train_df, test_df, model, model_dict, scaler = inputoutput_files(mod_name, path_model)


            train_df_smiles = train_df['SMILES'] 
            y_train_obs = train_df['observed'] 
    
            test_df_smiles = test_df['SMILES'] 
            y_test_obs = test_df['observed'] 
    
            print(f'\t\t Model type: {model._estimator_type}')
            
            print('\t[++] Calculating original metrics')
            
            old_new = 'old'
            
            df_metrics_byone = getMetrics(mod_name, model, train_df, test_df, old_new)
            
            df_metrics = pd.concat([df_metrics,df_metrics_byone], axis = 0)

            print('\t[++] Refitting model')
    
            
            
            
            if calculate == 'WITH':
            
                
                print('\t\t[+++] Recalculating descriptors')
        
                
                train_with_descriptors = calc_descriptors(option, mod_name, train_df, model_dict, 'train', results_folder)
                test_with_descriptors = calc_descriptors(option, mod_name, test_df, model_dict, 'test', results_folder)
                
            else:
                
                print('\t\t[+++] Loading new calculated descriptors')
                
                train_with_descriptors = pd.read_csv(results_folder + os.path.sep + f'raw_{mod_name}-descriptors-train.txt', sep ='\t')
                test_with_descriptors = pd.read_csv(results_folder + os.path.sep + f'raw_{mod_name}-descriptors-test.txt', sep ='\t')
                
    
    
            print('\t\t[+++] Imputing')
            
            if 'SMILES' in train_with_descriptors.columns:
                train_with_descriptors = train_with_descriptors.iloc[:, 1:]
                test_with_descriptors = test_with_descriptors.iloc[:, 1:]
    
            set_config(working_memory=4)
        
            train_descriptors_nanyfied = nanify(train_with_descriptors)
    
            imputer_estimator = KNNImputer(missing_values=np.nan,
                                 n_neighbors=3, weights="uniform")
    
            imputer = imputer_estimator.fit(train_with_descriptors)
            set_config(working_memory=1024)
    
    
            test_descriptors_nanyfied = nanify(test_with_descriptors)
            
            imputed_train_df = impute(train_descriptors_nanyfied, imputer)
            imputed_test_df = impute(test_descriptors_nanyfied, imputer)
            
            ##@@
            imputed_test_as_df = pd.DataFrame(imputed_test_df)
            imputed_test_as_df.columns = test_descriptors_nanyfied.columns
            imputed_test_as_df.to_csv(results_folder + os.path.sep + f'EX_{mod_name}-test_imputed.txt', sep = '\t', index = False)
            ##@@

            
            imputer_file_name = results_folder + os.path.sep + f'{mod_name}-train_set.imp'
            pickle.dump(imputer, open(imputer_file_name, 'wb'))
            print(f'\t\t\t imputer file created: {imputer_file_name}')
            
            
            print('\t\t[+++] Scaling')
            
            if isinstance(scaler, StandardScaler):
                print(f"\t\t\tusing {mod_name}-train_set.sca as StandardScaler")
                new_scaler = StandardScaler()
            elif isinstance(scaler, MinMaxScaler):
                print(f"\t\t\tusing {mod_name}-train_set.sca as MinMaxScaler")
                new_scaler = MinMaxScaler()
            else:
                print("Scaler file not found.")
                print("StandardScaler used by default")
                new_scaler = StandardScaler()        
            
            fitted_scaler = new_scaler.fit(imputed_train_df)
            
            scaler_file_name = results_folder + os.path.sep + f'{mod_name}-train_set.sca'
            pickle.dump(fitted_scaler, open(scaler_file_name, 'wb'))
            print(f'\t\t\t scaler file created: {scaler_file_name}')
            
            imputed_scaled_train = scale(imputed_train_df, fitted_scaler)
            imputed_scaled_train_df = pd.DataFrame(imputed_scaled_train)
            imputed_scaled_train_df.columns = train_with_descriptors.columns
            
            
            imputed_scaled_test = scale(imputed_test_df, fitted_scaler)
            imputed_scaled_test_df  = pd.DataFrame(imputed_scaled_test)
            imputed_scaled_test_df.columns = test_with_descriptors.columns
            
            ##@@
            imputed_scaled_test_df.to_csv(results_folder + os.path.sep + f'EX_{mod_name}-test_scaled.txt', sep = '\t', index = False)
            ##@@
            print('\t\t[+++] Fitting new model and predicting')
            
            fitted_model = model.fit(imputed_scaled_train_df, y_train_obs)
            
            model_file_name = results_folder + os.path.sep + f'{mod_name}.sav'
            pickle.dump(fitted_model, open(model_file_name, 'wb'))
            print(f'\t\t\t new model file created: {model_file_name}')        
            
            
            y_pred_train = predict(imputed_scaled_train_df, fitted_model)
            y_pred_test = predict(imputed_scaled_test_df, fitted_model)
            
            imputed_scaled_train_df.insert(0, 'SMILES', train_df_smiles)
            imputed_scaled_train_df.insert(len(imputed_scaled_train_df.columns), 'observed', y_train_obs)
            imputed_scaled_train_df.insert(len(imputed_scaled_train_df.columns), 'predicted', y_pred_train)
            
            imputed_scaled_test_df.insert(0, 'SMILES', test_df_smiles)
            imputed_scaled_test_df.insert(len(imputed_scaled_test_df.columns), 'observed', y_test_obs)
            imputed_scaled_test_df.insert(len(imputed_scaled_test_df.columns), 'predicted', y_pred_test)
    
    
            
            descriptors_file_name_train = results_folder + os.path.sep + f'{mod_name}-descriptors-train.txt'
            imputed_scaled_train_df.to_csv(descriptors_file_name_train, sep = '\t', index = False)
            descriptors_file_name_test = results_folder + os.path.sep + f'{mod_name}-descriptors-test.txt'
            imputed_scaled_test_df.to_csv(descriptors_file_name_test, sep = '\t', index = False)
            print(f'\t\t\t new file with descriptors created: {descriptors_file_name_train}')
            print(f'\t\t\t new file with descriptors created: {descriptors_file_name_test}')
            
            json_file_name = results_folder + os.path.sep + f'{mod_name}.json'
            with open(json_file_name, 'w', encoding='utf8') as json_file:
                json.dump(model_dict, json_file, ensure_ascii=False, indent=4)
            print(f'\t\t\t new json file created: {json_file_name}')
                
            print('\t[++] Calculating new metrics')
            
            old_new = 'new'
            
            df_metrics_byone = getMetrics(mod_name, model, imputed_scaled_train_df, imputed_scaled_test_df, old_new)
            
            individual_output_metrics_filename =   f'metrics_{mod_name}_{time.strftime("%Y%m%d_%H%M%S")}.csv'  
            print(f'\t\t metrics file created: {individual_output_metrics_filename}')
            df_metrics_byone.to_csv(individual_output_metrics_filename, sep = ';') 
            
            df_metrics = pd.concat([df_metrics,df_metrics_byone], axis = 0)    
                
                
                
                
        output_metrics_filename =   f'metrics_{time.strftime("%Y%m%d_%H%M%S")}.csv'  
        print(f'\t\t metrics file created: {output_metrics_filename}')
        df_metrics.to_csv(output_metrics_filename, sep = ';')    
                

    elif option == 'calculate':
        
        for mod_name in df_models_to_process['final_model_name']:
        
            model_data_folder =  '..' + os.path.sep + '..' + os.path.sep + 'Models'
            path_model = model_data_folder + os.path.sep + mod_name + os.path.sep + f'imputed-{mod_name}'
            
            topredict_folder = '.' + os.path.sep + 'predictions'
            
            model_dict =  load_json(path_model + os.path.sep + f'{mod_name}.json')
            
            for row, value in df_df_to_process.iterrows():
                
                dataset_name = value['dataset_iD']
                
                print(f'\t\t[+++] Working on {dataset_name} dataframe')
                
                df_to_calculate = pd.read_csv(topredict_folder + os.path.sep +  f'{value["file_to_calculate"]}' , sep = value['sep_for_calculation'])
            
                df_calculated = calc_descriptors(option, mod_name, df_to_calculate, model_dict, dataset_name, topredict_folder)


            
    elif option == 'predict':
        
        for mod_name in df_models_to_process['final_model_name']:
            
            print(f'\n[+] Working on {mod_name} model')
            
            print('\t[++] Reading input model files')
 

            model_data_folder =  '..' + os.path.sep + '..' + os.path.sep + 'Models' 
            path_model = model_data_folder + os.path.sep + mod_name + os.path.sep + f'imputed-{mod_name}'
            
            topredict_folder = '.' + os.path.sep + 'predictions' 
            predictions_folder = '.' + os.path.sep + 'predictions' + os.path.sep + 'df_predicted'
    
            makedir(predictions_folder)
               
            train_df, test_df, model, model_dict, scaler = inputoutput_files(mod_name, path_model)
            
                        
            model_type = model._estimator_type
            
            imputer = pickle.load(open(path_model + os.path.sep + f'{mod_name}-train_set.imp', "rb"))

            for row, value in df_df_to_process.iterrows():
                
                
                dataset_name = value['dataset_iD']
                
                print(f'\t\t[+++] Working on {dataset_name} dataframe')

                
                df_to_predict = pd.read_csv(topredict_folder + os.path.sep +  f'raw_{mod_name}-descriptors-{dataset_name}.txt', sep = '\t')
                # df_to_predict = pd.read_csv(path_model + os.path.sep +  f'raw_{mod_name}-descriptors-train.txt', sep = '\t')
               
                descriptors_model = [descriptor.split('$')[1] for descriptor in list(model_dict.keys())]
                
                dfpredict_smiles = df_to_predict['SMILES']
                
                descriptors_dfpredict = df_to_predict[descriptors_model]
                
                print('\t\t\t[++++] Imputing')
                
                descriptors_dfpredict_nanyfied = nanify(descriptors_dfpredict)
                imputed_dfpredict = impute(descriptors_dfpredict_nanyfied, imputer)

                ##@@
                imputed_dfpredict_as_df = pd.DataFrame(imputed_dfpredict)
                imputed_dfpredict_as_df.columns = descriptors_dfpredict_nanyfied.columns
                imputed_dfpredict_as_df.to_csv(topredict_folder + os.path.sep + f'EX_{mod_name}-new_imputed.txt', sep = '\t', index = False)
                ##@@

                
                print('\t\t\t[++++] Scaling')
                
                imputed_scaled_dfpredict = scale(imputed_dfpredict, scaler)
                imputed_scaled_dfpredict_df  = pd.DataFrame(imputed_scaled_dfpredict)
                imputed_scaled_dfpredict_df.columns = descriptors_dfpredict.columns

                ##@@

                imputed_scaled_dfpredict_df.to_csv(topredict_folder + os.path.sep + f'EX_{mod_name}-new_scaled.txt', sep = '\t', index = False)
                ##@@

                
                ##@@
                
                ##@@               
                print('\t\t\t[++++] Predicting')
                
                y_pred_dfpredict = predict(imputed_scaled_dfpredict_df, model)
    
                imputed_scaled_dfpredict_df.insert(0, 'SMILES', dfpredict_smiles)
                imputed_scaled_dfpredict_df.insert(len(imputed_scaled_dfpredict_df.columns), 'predicted', y_pred_dfpredict)
                
                predicted_filename = predictions_folder + os.path.sep + f'{dataset_name}-{mod_name}-onlypredicted.csv'
                imputed_scaled_dfpredict_df.to_csv(predicted_filename, sep = ';', index = False)
                print(f'\t\t\t new file with descriptors created: {predicted_filename}')
                
                print('\t\t\t[++++] Extracting experimental values')
                
                combined_train_test = pd.concat([train_df, test_df], axis = 0)
                
                combined_train_test.set_index('SMILES', inplace = True)
                
                combined_train_test_dict = combined_train_test.T.to_dict()
                
                experimentals = []
                
                for smi in imputed_scaled_dfpredict_df['SMILES']:
                    if smi in combined_train_test_dict.keys():
                        if model_type == 'classifier':
                            experimentals.append(int(combined_train_test_dict[smi]['observed']))
                        else:
                            experimentals.append(combined_train_test_dict[smi]['observed'])
                    else:
                        experimentals.append(np.nan)
                        
                imputed_scaled_dfpredict_df.insert(len(imputed_scaled_dfpredict_df.columns), 'experimental', experimentals)
                
                print('\t\t\t[++++] Calculating AD')
                
                df_predict_with_ads = calculate_ADs(train_df,imputed_scaled_dfpredict_df,descriptors_model, mod_name)
                
                
                final_predicted_filename = predictions_folder + os.path.sep + f'{dataset_name}-{mod_name}-predictedAD.csv'
                df_predict_with_ads.to_csv(final_predicted_filename, sep = ';', index = False)
                print(f'\t\t\t new file with descriptors created: {final_predicted_filename}')

        
    elif option == 'mergeJSONS':
        
        merged_dict = {}
        
        counter = 0
        
        for mod_name in df_models_to_process['final_model_name']:
            
            model_data_folder =  '..' + os.path.sep + '..' + os.path.sep + 'Models' 
            path_model = model_data_folder + os.path.sep + mod_name + os.path.sep + f'imputed-{mod_name}'
            
            model_dict =  load_json(path_model + os.path.sep + f'{mod_name}.json')
            
            merged_dict.update(model_dict)
            
            counter = counter + len(model_dict)
            
        print(counter, len(merged_dict))
            
            
            

#%% ADs
# PATH2 = r'C:\Users\proto\Desktop\tothexinxol\IRB\Models\TK_OATP1B1inh_no3D_us\predictions'    
    
  
# train_df = pd.read_csv(f"{PATH2}/{mod_name}-descriptors-train.txt",sep="\t")    
# test_df = pd.read_csv(f"{PATH2}/predict_{mod_name}.csv",sep=";")  
# ads = calculate_ADs(train_df,test_df,mod_name)

# ads.to_csv(f"{PATH2}/predict_{mod_name}_ad.csv",sep=";",index=False)

'''
elapsed_time = 11.042634963989258 30 mols
'''

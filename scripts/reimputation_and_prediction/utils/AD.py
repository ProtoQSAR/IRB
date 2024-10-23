# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 16:16:33 2022

@author: Rita
"""

#script to get AD evaluated by Pravin's method and by leverage for screenings
#inputs: datasets with all descriptors used to train data
#        datasets resulting of screening


import pandas as pd

import numpy as np
import statistics as stats


from rdkit import Chem, DataStructs
from sklearn.metrics.pairwise import euclidean_distances
from rdkit.Chem import rdMolDescriptors





def pravin_AD_calculation_test(train_descriptors, test_descriptors, descriptors):

    #get mean and sd of train data for each descriptor and save them in two lists
    means_descriptors = []
    sd_descriptors = []
    for descriptor in descriptors:
        train_desc_x = train_descriptors.loc[:,descriptor]
        #compute mean in train dataset
        mean_desc = stats.mean(train_desc_x)
        means_descriptors.append(mean_desc)
        #compute standard deviation in train dataset
        sd_desc = stats.stdev(train_desc_x)
        sd_descriptors.append(sd_desc)

    #now inside each predictions file (test file), iterate by compound (by row)
    AD_results = []
    for index, compound in test_descriptors.iterrows():
        standarized_descriptors_compound = []
        for i, descriptor in enumerate(descriptors):
            #get standarized descriptor for each compound and each descriptor
            xki = compound[i]
            xi = means_descriptors[i]
            sdi = sd_descriptors[i]
            ski = abs(xki-xi)/sdi
            standarized_descriptors_compound.append(ski)
        if max(standarized_descriptors_compound) <= 3:
            AD = "In"
        elif min(standarized_descriptors_compound) >3:
            AD = "Out"
        else:
            skmean = stats.mean(standarized_descriptors_compound)
            sksd = stats.stdev(standarized_descriptors_compound)
            snew = skmean + (1.28*sksd)
            if snew > 3:
                AD = "Out"
            else:
                AD = "In"
        AD_results.append(AD)

    test_AD = AD_results

    #now inside each predictions file (test file), iterate by compound (by row)
    AD_results = []
    for index, compound in train_descriptors.iterrows():
        standarized_descriptors_compound = []
        for i, descriptor in enumerate(descriptors):
            #get standarized descriptor for each compound and each descriptor
            xki = compound[i]
            xi = means_descriptors[i]
            sdi = sd_descriptors[i]
            ski = abs(xki-xi)/sdi
            standarized_descriptors_compound.append(ski)
        if max(standarized_descriptors_compound) <= 3:
            AD = "In"
        elif min(standarized_descriptors_compound) >3:
            AD = "Out"
        else:
            skmean = stats.mean(standarized_descriptors_compound)
            sksd = stats.stdev(standarized_descriptors_compound)
            snew = skmean + (1.28*sksd)
            if snew > 3:
                AD = "Out"
            else:
                AD = "In"
        AD_results.append(AD)

    train_AD = AD_results

    return test_AD, train_AD

#function to get the leverage value for each compound in an array
def get_leverage(descriptor_array, ainv):
    leverages = []
    for row in descriptor_array:
        xi = row
        xit = xi.transpose()
        lev = xi.dot(ainv).dot(xit)
        leverages.append(lev)

    return leverages

#function to transform leverages to "in" and "out"
def get_leverage_tags(leverages, hlim):
    tags = []
    for lev in leverages:
        if lev > hlim:
            tags.append("Out")
        else:
            tags.append("In")

    return tags

#from ProtoPRED get fingerprints
def get_fingerprints_rdkit(mol): # NEW DESCRIPTORS USED
    '''Auxiliary method for calculating rdkit fingerprints.
    '''
    try:
        # print(Chem.MolToSmiles(mol), rdMolDescriptors.GetMACCSKeysFingerprint(mol))
        return rdMolDescriptors.GetMACCSKeysFingerprint(mol)

    except:
        pass

#from ProtoPRED tanimoto AD
def tanimoto(train_fps, smiles, top_scores=10):
    '''Calculates Tanimoto similarity between query SMILES and training set.
    '''
    results = []
    result_tanimoto_dict = {}
    # complete_dict = {smiles:{}}

    try:
        # print(smiles)
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        fp = get_fingerprints_rdkit(mol)
        # print('[+++]  Calculating Tanimoto coefficient')
        # print(self.train_fps)
        for smi,fp_train in train_fps.items():
            try:
                res = round(DataStructs.FingerprintSimilarity(fp, fp_train), 2)
            # Morgan fingerprints not correctly calculated
            except:
                res = 0.0

            # print("tanimoto: ", res)

            result_tanimoto_dict[smi] = res
            results.append(res)

        results.sort(reverse=True)
    except:
        results.append(0.0)

    return results[:top_scores],result_tanimoto_dict

#from ProtoPRED euclidean distance AD
def euclidean_distance(X_train, smis_train, X_test, smis_test):
    '''Estimates the AD based on the euclidean distance between the training
    set descriptors and the query molecule
    '''

    ed_criteria = []

    # calculate ed distance within train dataset
    ed_all = euclidean_distances(X_train, X_train)

    ed_all.sort(axis=1)
    # print('ed_all',ed_all)
    last = ed_all[:,-1]
    # print('last',last)
    threshold_ed_train = max(last)
    # print('threshold_ed_train',threshold_ed_train)# threshold for max distance

    dictionary_euclidean = {}
    euclidean_complete_dict = {}
    euclidean_query = []
    # print(self.train_desc.iloc[0])
    # print(self.query_descriptors)
    # calculate ed distance with query mols and train dataset

    for query_smile, row in zip(smis_test,X_test):
        # print(row)
        row = row.reshape(1, row.shape[0])
        ed =  euclidean_distances(X_train, row)

        edlisst = list(ed.flatten())

        for smiles, distance in zip(smis_train,edlisst):
            dictionary_euclidean[smiles] = distance.round(4)

        euclidean_complete_dict[query_smile] = dictionary_euclidean

        max_ed_molecule = max(edlisst)

        # print('max_ed_molecule',max_ed_molecule)
        # print('dictionary_euclidean',dictionary_euclidean)

        euclidean_query.append(max_ed_molecule.round(4))

        if max_ed_molecule < threshold_ed_train:
            ed_criteria.append("In")
        else:
            ed_criteria.append("Out")

    return ed_criteria, euclidean_complete_dict,euclidean_query

# from ProtoPRED range AD
def get_range_ad(df_descr_train, df_descr_test):

    train_maxdesc = []
    train_mindesc = []
    for col in df_descr_train.columns:
        train_maxdesc.append(df_descr_train[col].max())
        train_mindesc.append(df_descr_train[col].min())


    range_ad = []
    range_bad = []
    for index, row in df_descr_test.iterrows():
        bad=[]
        for i, col in enumerate(df_descr_test.columns):
            if row[col] > train_maxdesc[i] or row[col] < train_mindesc[i]:
                bad.append({"name":col,
                            "value":row[col],
                            "min":train_mindesc[i],
                            "max":train_maxdesc[i]})
        range_bad.append(bad)

        if len(bad) > 0:
            range_ad.append("Out")
        else:
            range_ad.append("In")

    return range_ad, range_bad



def calculate_ADs(train_df: pd.DataFrame,
                    test_df: pd.DataFrame,
                    descriptors_model: list,
                    model = str
                    ):
    train_fps = {}
    for smiles_train in train_df['SMILES'].values:
        # print(smiles_train)
        mol = Chem.MolFromSmiles(smiles_train, sanitize=True)
        fp = get_fingerprints_rdkit(mol)
        fp
        train_fps[smiles_train] = fp

    #get dfs of descriptors for train and test set and a list of descriptors for the model
    test_descriptors = test_df[descriptors_model]
    
    train_descriptors = train_df[descriptors_model]

    #get two lists with "in" and "out" according to Pravin's AD calculation method
    # test_AD, train_AD = pravin_AD_calculation_test(train_descriptors, test_descriptors, descriptors)

    ###let's calculate leverage for test and train for this model
    #convert df to numpy
    train_descriptor_array = train_descriptors.to_numpy()
    test_descriptor_array = test_descriptors.to_numpy()
    #make the matrix calculations to obtain A = (XtX)^-1
    xt_array_train = train_descriptor_array.transpose()
    a_array_train = xt_array_train.dot(train_descriptor_array)
    a_inv_train = np.linalg.inv(a_array_train)

    #save train leverages in a list
    train_leverages = get_leverage(train_descriptor_array, a_inv_train)
    #save test leverages in a list
    test_leverages =  get_leverage(test_descriptor_array, a_inv_train)

    #get h*
    n = len(train_leverages)
    p = len(train_descriptors.columns)
    h_limit = 3*p/n

    #convert leverage values to in/out
    # train_lev_tags = get_leverage_tags(train_leverages, h_limit)
    test_lev_tags = get_leverage_tags(test_leverages, h_limit)

    #get y prediction from test file
    # test_y_predict = test_df.iloc[:,-1]

    test_y_predict = test_df['predicted']
    test_y_exp = test_df['experimental']

    #get TANIMOTO from test
    tanimoto_best_scores = []
    tanimoto_complete_dict = {}

    for smi in test_df['SMILES'].values:

        best_tanimotos,dict_tanimotos =  tanimoto(train_fps, smi, top_scores=1)
        tanimoto_complete_dict[smi] = dict_tanimotos
        tanimoto_best_scores.append(best_tanimotos)

    tanimotos = [tupla[0] for tupla in tanimoto_best_scores]
    threshold = 0.528 # this is the Tanimoto similarity threshold for MACCS for a 90% level (https://rdkit.blogspot.com/2013/10/fingerprint-thresholds.html)

    tanimoto_inout = ["In" if tanimotos[i] > threshold else 'Out'
                for i in range(len(tanimotos))]

    #get EUCLIDEAN DISTANCE AD

    euclidean_distance_ad, euclidean_complete_dict,eucledian_query = euclidean_distance(train_descriptor_array,
                                                                                        train_df['SMILES'].values,
                                                                                        test_descriptor_array,
                                                                                        test_df['SMILES'].values)

    #get RANGE AD
    range_ad, range_bad = get_range_ad(train_descriptors, test_descriptors)

    #generate a dataframe with data from test
    model_results_df = pd.DataFrame({"{} predicted value".format(model) : test_y_predict,
                                        "{} experimental value".format(model): test_y_exp,
                                        # "{} AD leverage2".format(model) : test_AD,
                                        "{} AD leverage".format(model) : test_lev_tags,
                                        "{} AD tanimoto".format(model) : tanimoto_inout,
                                        "{} AD euclidean".format(model) : euclidean_distance_ad,
                                        "{} AD range".format(model) : range_ad,
                                        })

    conditions = [
        model_results_df["{} AD tanimoto".format(model)].eq("In") &
        model_results_df["{} AD euclidean".format(model)].eq("In") &
        model_results_df["{} AD leverage".format(model)].eq("In") &
        # model_results_df["{} AD leverage2".format(model)].eq("In") &
        model_results_df["{} AD range".format(model)].eq("In"),
        # model_results_df["{} AD leverage2".format(model)].eq("In"),
        model_results_df["{} AD leverage".format(model)].eq("In"),
        model_results_df["{} AD tanimoto".format(model)].eq("In"),
        model_results_df["{} AD euclidean".format(model)].eq("In"),
        model_results_df["{} AD range".format(model)].eq("In"),

        ]

    choices =  [
        "All in",
        # "Leverage2",
        "Leverage",
        "Tanimoto",
        "Euclidean",
        "Range"
        ]

    model_results_df["AD"] = np.select(conditions, choices, default="Out")

    # obtain AD column as in ProtoPRED

    dict_protopred_ad_tags = {f"{model} AD tanimoto":"T",
                                f"{model} AD euclidean":"E",
                                f"{model} AD leverage":"L",
                                # f"{model} AD leverage2":"L2",
                                f"{model} AD range":"R"}

    # function to get which AD columns are "In" for each compound
    def get_ad_tags(row):
        ad_tags = []
        for key, value in dict_protopred_ad_tags.items():
            if row[key] == "In":
                ad_tags.append(value)

        # generate string for total AD column
        if len(ad_tags) == 0:
            ad_col_tag = "Outside"
        else:
            # generate string with AD tags
            ad_tags = "/".join(ad_tags)
            ad_col_tag = f"Inside ({ad_tags})"
        return ad_col_tag

    model_results_df["AD_ProtoPRED"] = model_results_df.apply(get_ad_tags, axis=1)




    test_sub = test_df.iloc[:,:-2]
    results_join = pd.concat([test_sub,model_results_df], axis=1)

    return results_join


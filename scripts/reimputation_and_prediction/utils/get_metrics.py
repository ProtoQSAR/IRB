# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 12:51:21 2022

@author: EvaSC
"""
import os

import pickle
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.metrics import *
import numpy as np
import pandas as pd

def classifier_assessment(Y_train, Y_train_pred):
    # print(Y_train,Y_train_pred)
    # WRITE YOUR CODE HERE
    data_length = len(Y_train)

    FP = 0
    FN = 0
    TP = 0
    TN = 0
    IP = 0
    IN = 0

    for i in range(data_length):
        if Y_train[i] == 0 and Y_train_pred[i] == 0:
            TN = TN+1
        elif Y_train[i] == 1 and Y_train_pred[i] == 0:
            FN = FN+1
        elif Y_train[i] == 1 and Y_train_pred[i] == 1:
            TP = TP+1
        elif Y_train[i] == 0 and Y_train_pred[i] == 1:
            FP = FP+1
        elif Y_train[i] == 0 and Y_train_pred[i] == 2:
            IN = IN+1
        elif Y_train[i] == 1 and Y_train_pred[i] == 2:
            IP = IP+1

    try:
        accuracy = (TP + TN) / (TP + TN + FP + FN)
    except:
        accuracy = -1.0
    try:
        sensitivity = TP / (TP + FN)
    except:
        sensitivity = -1.0
    try:
        specificity = TN / (TN + FP)
    except:
        specificity = -1.0
    try:
        precision = TP / (TP + FP)
    except:
        precision= -1.0

    accuracy, sensitivity, specificity, precision = round(accuracy, 2), round(
        sensitivity, 2), round(specificity, 2), round(precision, 2)

    return [accuracy, sensitivity, specificity, precision]

def classifier_assessment_cv(together_for_observations, model):



    X_together = together_for_observations.iloc[:, 1:-2]

    y_together = together_for_observations['observed']


    def confusion_matrix_metrics(clf, X, y):
        y_pred = clf.predict(X)
        cm = confusion_matrix(y, y_pred)
        tn, fn, tp, fp = cm[0, 0], cm[1, 0], cm[1, 1], cm[0, 1]
        accuracy = (tn+tp)/(tn+tp+fn+fp)
        precision = tp/(tp+fp)
        recall = tp/(tp+fn)
        f1_score = 2*precision*recall/(precision+recall)
        sensit = tp/(tp+fn)
        spec = tn/(tn+fp)
        NPV = tn/(tn+fn)
        FNR = fn/(fn+tp)
        FPR = fp/(fp+tn)
        FDR = fp/(fp+tp)
        FOR = fn/(fn+tn)
        F_score = 2*tp/(2*tp+fp+fn)
        MCC = ((tp*tn)-(fp*fn))/(np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
        CSI = tp/(tp+fn+fp)

        return {'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1_score': f1_score,
                'sensit': sensit, 'spec': spec, 'NPV': NPV, 'FNR': FNR, 'FPR': FPR, 'FDR': FDR,
                'FOR': FOR, 'F_score': F_score, 'MCC': MCC, 'CSI': CSI}

    scores_together = cross_validate(model, X_together, y_together, cv=5,
                            scoring=confusion_matrix_metrics, return_train_score=True)


    scores_together_auc = cross_validate(model, X_together, y_together, cv=5,
                            scoring='roc_auc', return_train_score=True)

    CV_acc_train_mean = np.mean(scores_together['train_accuracy'])
    CV_prec_train_mean = np.mean(scores_together['train_precision'])
    CV_recall_train_mean = np.mean(scores_together['train_recall'])
    CV_f1_train_mean = np.mean(scores_together['train_f1_score'])
    CV_roc_auc_train_mean = np.mean(scores_together_auc['train_score'])
    CV_sensit_train_mean = np.mean(scores_together['train_sensit'])
    CV_spec_train_mean = np.mean(scores_together['train_spec'])
    CV_NPV_train_mean = np.mean(scores_together['train_NPV'])
    CV_FNR_train_mean = np.mean(scores_together['train_FNR'])
    CV_FPR_train_mean = np.mean(scores_together['train_FPR'])
    CV_FDR_train_mean = np.mean(scores_together['train_FDR'])
    CV_FOR_train_mean = np.mean(scores_together['train_FOR'])
    CV_F_score_train_mean = np.mean(scores_together['train_F_score'])
    CV_MCC_train_mean = np.mean(scores_together['train_MCC'])
    CV_CSI_train_mean = np.mean(scores_together['train_CSI'])

    CV_acc_train_SD = np.std(scores_together['train_accuracy'])
    CV_prec_train_SD = np.std(scores_together['train_precision'])
    CV_recall_train_SD = np.std(scores_together['train_recall'])
    CV_f1_train_SD = np.std(scores_together['train_f1_score'])
    CV_roc_auc_train_SD = np.std(scores_together_auc['train_score'])
    CV_sensit_train_SD = np.std(scores_together['train_sensit'])
    CV_spec_train_SD = np.std(scores_together['train_spec'])
    CV_NPV_train_SD = np.std(scores_together['train_NPV'])
    CV_FNR_train_SD = np.std(scores_together['train_FNR'])
    CV_FPR_train_SD = np.std(scores_together['train_FPR'])
    CV_FDR_train_SD = np.std(scores_together['train_FDR'])
    CV_FOR_train_SD = np.std(scores_together['train_FOR'])
    CV_F_score_train_SD = np.std(scores_together['train_F_score'])
    CV_MCC_train_SD = np.std(scores_together['train_MCC'])
    CV_CSI_train_SD = np.std(scores_together['train_CSI'])

    CV_acc_test_mean = np.mean(scores_together['test_accuracy'])
    CV_prec_test_mean = np.mean(scores_together['test_precision'])
    CV_recall_test_mean = np.mean(scores_together['test_recall'])
    CV_f1_test_mean = np.mean(scores_together['test_f1_score'])
    CV_roc_auc_test_mean = np.mean(scores_together_auc['test_score'])
    CV_sensit_test_mean = np.mean(scores_together['test_sensit'])
    CV_spec_test_mean = np.mean(scores_together['test_spec'])
    CV_NPV_test_mean = np.mean(scores_together['test_NPV'])
    CV_FNR_test_mean = np.mean(scores_together['test_FNR'])
    CV_FPR_test_mean = np.mean(scores_together['test_FPR'])
    CV_FDR_test_mean = np.mean(scores_together['test_FDR'])
    CV_FOR_test_mean = np.mean(scores_together['test_FOR'])
    CV_F_score_test_mean = np.mean(scores_together['test_F_score'])
    CV_MCC_test_mean = np.mean(scores_together['test_MCC'])
    CV_CSI_test_mean = np.mean(scores_together['test_CSI'])

    CV_acc_test_SD = np.std(scores_together['test_accuracy'])
    CV_prec_test_SD = np.std(scores_together['test_precision'])
    CV_recall_test_SD = np.std(scores_together['test_recall'])
    CV_f1_test_SD = np.std(scores_together['test_f1_score'])
    CV_roc_auc_test_SD = np.std(scores_together_auc['test_score'])
    CV_sensit_test_SD = np.std(scores_together['test_sensit'])
    CV_spec_test_SD = np.std(scores_together['test_spec'])
    CV_NPV_test_SD = np.std(scores_together['test_NPV'])
    CV_FNR_test_SD = np.std(scores_together['test_FNR'])
    CV_FPR_test_SD = np.std(scores_together['test_FPR'])
    CV_FDR_test_SD = np.std(scores_together['test_FDR'])
    CV_FOR_test_SD = np.std(scores_together['test_FOR'])
    CV_F_score_test_SD = np.std(scores_together['test_F_score'])
    CV_MCC_test_SD = np.std(scores_together['test_MCC'])
    CV_CSI_test_SD = np.std(scores_together['test_CSI'])

    print('\n\t\t######################################################')
    print("\t\t############## Cross-validation results ##############")
    print('\t\t######################################################')

    print('\n\t\t\t|Train          | Test')
    print('\t\t    acc\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_acc_train_mean,CV_acc_train_SD,CV_acc_test_mean,CV_acc_test_SD))
    print('\t\t   prec\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_prec_train_mean,CV_prec_train_SD,CV_prec_test_mean,CV_prec_test_SD))
    print('\t\t recall\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_recall_train_mean,CV_recall_train_SD,CV_recall_test_mean,CV_recall_test_SD))
    print('\t\t     f1\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_f1_train_mean,CV_f1_train_SD,CV_f1_test_mean,CV_f1_test_SD))
    print('\t\t    auc\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_roc_auc_train_mean,CV_roc_auc_train_SD,CV_roc_auc_test_mean,CV_roc_auc_test_SD))
    print('\t\t sensit\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_sensit_train_mean,CV_sensit_train_SD,CV_sensit_test_mean,CV_sensit_test_SD))
    print('\t\t   spec\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_spec_train_mean,CV_spec_train_SD,CV_spec_test_mean,CV_spec_test_SD))
    print('\t\t    NPV\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_NPV_train_mean,CV_NPV_train_SD,CV_NPV_test_mean,CV_NPV_test_SD))
    print('\t\t    FNR\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_FNR_train_mean,CV_FNR_train_SD,CV_FNR_test_mean,CV_FNR_test_SD))
    print('\t\t    FPR\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_FPR_train_mean,CV_FPR_train_SD,CV_FPR_test_mean,CV_FPR_test_SD))
    print('\t\t    FDR\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_FDR_train_mean,CV_FDR_train_SD,CV_FDR_test_mean,CV_FDR_test_SD))
    print('\t\t    FOR\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_FOR_train_mean,CV_FOR_train_SD,CV_FOR_test_mean,CV_FOR_test_SD))
    print('\t\tF_score\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_F_score_train_mean,CV_F_score_train_SD,CV_F_score_test_mean,CV_F_score_test_SD))
    print('\t\t    MCC\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_MCC_train_mean,CV_MCC_train_SD,CV_MCC_test_mean,CV_MCC_test_SD))
    print('\t\t    CSI\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_CSI_train_mean,CV_CSI_train_SD,CV_CSI_test_mean,CV_CSI_test_SD))


def regressor_assessment(Y_train, Y_train_pred):

    rsq = r2_score(Y_train, Y_train_pred)
    ev = explained_variance_score(Y_train, Y_train_pred)
    mse = mean_squared_error(Y_train, Y_train_pred)
    mae = mean_absolute_error(Y_train, Y_train_pred)

    # Keep the prints to check your work
    rsq, ev, mse, mae = round(rsq, 2), round(
        ev, 2), round(mse, 2), round(mae, 2)

    return [rsq, ev, mse, mae]

def regressor_assessment_cv(together_for_observations, model):

    X_together = together_for_observations.iloc[:, 1:-2]


    y_together = together_for_observations['observed']


    scores_together = cross_validate(model, X_together, y_together, cv=5,
                            scoring=( 'explained_variance','neg_mean_absolute_error',
                            'neg_mean_squared_error','neg_median_absolute_error','r2'),
                            return_train_score=True)


    CV_evs_train_mean = np.mean(scores_together['train_explained_variance'])
    CV_mae_train_mean = -np.mean(scores_together['train_neg_mean_absolute_error'])
    CV_mse_train_mean = -np.mean(scores_together['train_neg_mean_squared_error'])
    CV_medae_train_mean = -np.mean(scores_together['train_neg_median_absolute_error'])
    CV_r2_train_mean = np.mean(scores_together['train_r2'])

    CV_evs_train_SD = np.std(scores_together['train_explained_variance'])
    CV_mae_train_SD = np.std(scores_together['train_neg_mean_absolute_error'])
    CV_mse_train_SD = np.std(scores_together['train_neg_mean_squared_error'])
    CV_medae_train_SD = np.std(scores_together['train_neg_median_absolute_error'])
    CV_r2_train_SD = np.std(scores_together['train_r2'])


    CV_evs_test_mean = np.mean(scores_together['test_explained_variance'])
    CV_mae_test_mean = -np.mean(scores_together['test_neg_mean_absolute_error'])
    CV_mse_test_mean = -np.mean(scores_together['test_neg_mean_squared_error'])
    CV_medae_test_mean = -np.mean(scores_together['test_neg_median_absolute_error'])
    CV_r2_test_mean = np.mean(scores_together['test_r2'])

    CV_evs_test_SD = np.std(scores_together['test_explained_variance'])
    CV_mae_test_SD = np.std(scores_together['test_neg_mean_absolute_error'])
    CV_mse_test_SD = np.std(scores_together['test_neg_mean_squared_error'])
    CV_medae_test_SD = np.std(scores_together['test_neg_median_absolute_error'])
    CV_r2_test_SD = np.std(scores_together['test_r2'])

    print('\n\t\t######################################################')
    print("\t\t############## Cross-validation results ##############")
    print('\t\t######################################################')



    print('\n\t\t\t|Train          | Test')
    print('\t\t    EV\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_evs_train_mean,CV_evs_train_SD,CV_evs_test_mean,CV_evs_test_SD))
    print('\t\t   MAE\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_mae_train_mean,CV_mae_train_SD,CV_mae_test_mean,CV_mae_test_SD))
    print('\t\t   MSE\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_mse_train_mean,CV_mse_train_SD,CV_mse_test_mean,CV_mse_test_SD))
    print('\t\t MEDAE\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_medae_train_mean,CV_medae_train_SD,CV_medae_test_mean,CV_medae_test_SD))
    print('\t\tRscore\t|%3.2f +/- %3.3f |%3.2f  +/- %3.3f' %(CV_r2_train_mean,CV_r2_train_SD,CV_r2_test_mean,CV_r2_test_SD))




def getMetrics(mod_name, model, df_train, df_test, old_new):
    dictio_metrics = {}    
    model_type = model._estimator_type
    
    model_name = mod_name + f'_{old_new}'
    
    y_train = df_train['observed']
    y_train_pred = df_train['predicted']	

    y_test = df_test['observed']
    y_test_pred = df_test['predicted']

    merged = pd.concat([df_train,df_test])	  

    if model_type == 'classifier':
        train_metrics = classifier_assessment(y_train, y_train_pred)
        test_metrics = classifier_assessment(y_test, y_test_pred)
        classifier_assessment_cv(merged, model)
        
    elif model_type == 'regressor':
        train_metrics = regressor_assessment(y_train, y_train_pred)
        test_metrics = regressor_assessment(y_test, y_test_pred)
        regressor_assessment_cv(merged, model)        
        
    dictio_metrics[model_name] = {}
    dictio_metrics[model_name]['Model type'] = model_type

    dictio_metrics[model_name]['Train R2/Acc'] = train_metrics[0]
    dictio_metrics[model_name]['Train Ev/Sensit'] = train_metrics[1]
    dictio_metrics[model_name]['Train MSE/Specif'] = train_metrics[2]
    dictio_metrics[model_name]['Train MAE/Precis'] = train_metrics[3]

    dictio_metrics[model_name]['Test R2/Acc'] = test_metrics[0]
    dictio_metrics[model_name]['Test Ev/Sensit'] = test_metrics[1]
    dictio_metrics[model_name]['Test MSE/Specif'] = test_metrics[2]
    dictio_metrics[model_name]['Test MAE/Precis'] = test_metrics[3]

    df_metrics = pd.DataFrame.from_dict(dictio_metrics,  orient='index')  

    return df_metrics     
        

    # classifier_assessment_cv(together_for_observations, path, model)

    # if cv_only == False:

    #     train_metrics = regressor_assessment(y_train, y_train_pred2)
    #     test_metrics = regressor_assessment(y_test, y_test_pred2)
    # else:
    #     print('cv metrics')
    #     regressor_assessment_cv(together_for_observations, path, model)

if __name__ == '__main__':
    print('recalculating')
    df_recalc = getMetrics(model)

    # df_recalc = main(recalculate=False, cv_only=False)

    # print('nonrecalculating')
    # df_norecalc =main(recalculate=False)

    # print(df_recalc == df_norecalc)

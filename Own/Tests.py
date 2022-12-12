# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 11:17:12 2022

@author: JK-WORK
"""


import csv
import pandas as pd
import copy
import numpy as np


# df_rad = pd.read_csv('Dataset/radiometries_train.csv');
# df_rad_ids_ord=df_rad.copy()
# true_labels=df_rad["IDs"].to_numpy();
# label_values=np.unique(true_labels)
# label_values.sort()
# df_rad_ids_ord['IDs'] = df_rad_ids_ord['IDs'].apply(lambda x: np.argwhere(label_values==x)[0,0])
# true_labels_ord=df_rad["IDs"].to_numpy();
# unique, counts = np.unique(true_labels_ord, return_counts=True)
# prior_probs=counts/np.sum(counts);
# true_labels_norm = [np.argwhere(label_values==x)[0,0] for x in true_labels]

# df_rad_test = pd.read_csv('Dataset/radiometries_test.csv')
# df_rad_test_ids_ord=df_rad_test.copy()
# df_rad_test_ids_ord['IDs'] = df_rad_test_ids_ord['IDs'].apply(lambda x: np.argwhere(label_values==x)[0,0])
# true_labels_test=df_rad_test["IDs"].to_numpy();
# label_values_test=np.unique(true_labels_test)
# label_values_test.sort()

# pred_hybridHMM_test = pd.read_csv('Dataset/predictions_hybrid_test.csv')
# pred_hybridHMM_test['best_pep_iz'] = pred_hybridHMM_test['best_pep_iz'].apply(lambda x: np.argwhere(label_values==x)[0,0])
# HMM_Hybrid_y_pred=pred_hybridHMM_test['best_pep_iz'].to_numpy();

df_train = pd.read_csv('Dataset/radiometries_train.csv')
y_train=df_train["IDs"].to_numpy();
print("True-id array max value: "+str(np.max(y_train))+", and unique values: "+str(len(np.unique(y_train))))
df_train["IDs"].hist(bins=np.max(y_train))
df_train_ynorm=df_train.copy()
labels=np.unique(y_train)
labels.sort()
df_train_ynorm['IDs'] = df_train_ynorm['IDs'].apply(lambda x: np.argwhere(labels==x)[0,0])
# df_train_ynorm['IDs'].hist(bins=len(labels))

# df_train_ynorm_xint=df_train_ynorm.copy()
# df_train_ynorm_xint.iloc[:,:-1]=abs(np.rint((df_train.iloc[:,:-1].to_numpy()/10000)))
# df_train_ynorm_xint.iloc[:,:-1]=df_train_ynorm_xint.iloc[:,:-1].astype(int)


# n_edman=9;n_ch=3;
# first_meas_chs=[i*(n_edman+1) for i in range(n_ch)]

# temp=df_train_ynorm_xint.iloc[:,:-1].to_numpy();
# temp_channels=np.hsplit(temp, n_ch)
# temp_channels_diff=[np.concatenate((temp_channels[i][:,0].reshape(-1,1),np.diff(temp_channels[i],axis=1)),axis=1) for i in range(n_ch)]
# array_train_ynorm_xint_diff=copy.deepcopy(temp_channels_diff[0]);
# for i in range(1,n_ch):
#     array_train_ynorm_xint_diff=np.concatenate((array_train_ynorm_xint_diff,temp_channels_diff[i]),axis=1)
# df_train_ynorm_xint_diff=df_train_ynorm_xint.copy()
# df_train_ynorm_xint_diff.iloc[:,:-1]=array_train_ynorm_xint_diff;

df_test = pd.read_csv('Dataset/radiometries_test.csv')
y_test=df_test["IDs"].to_numpy();
y_test_norm = [np.argwhere(labels==x)[0,0] for x in y_test]

df_test_ynorm=df_test.copy()
df_test_ynorm['IDs'] = y_test_norm

def normalize(df):
    df_norm=df.copy()
    array_norm_wo_label=df_norm.iloc[:,:-1].to_numpy();
    mean=np.mean(array_norm_wo_label);
    std=np.std(array_norm_wo_label)
    norm=[mean, std]
    df_norm.iloc[:,:-1] = (df_norm.iloc[:,:-1] - mean)/std;
    return df_norm,norm;
def apply_norm(df, norm):
    mean,std=norm;
    df_norm=df.copy()
    df_norm.iloc[:,:-1] = (df_norm.iloc[:,:-1] - mean)/std;
    return df_norm;


df_train_ynorm_xnorm,norm=normalize(df_train_ynorm)
df_test_ynorm_xnorm=apply_norm(df_test_ynorm, norm)

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 11:17:12 2022

@author: JK-WORK
"""


import csv
import pandas as pd
import copy
import numpy as np


df_rad = pd.read_csv('Dataset/radiometries_train.csv');
df_rad_ids_ord=df_rad.copy()
true_labels=df_rad["IDs"].to_numpy();
label_values=np.unique(true_labels)
label_values.sort()
df_rad_ids_ord['IDs'] = df_rad_ids_ord['IDs'].apply(lambda x: np.argwhere(label_values==x)[0,0])
true_labels_ord=df_rad["IDs"].to_numpy();
unique, counts = np.unique(true_labels_ord, return_counts=True)
prior_probs=counts/np.sum(counts);


df_rad_test = pd.read_csv('Dataset/radiometries_test.csv')
df_rad_test_ids_ord=df_rad_test.copy()
df_rad_test_ids_ord['IDs'] = df_rad_test_ids_ord['IDs'].apply(lambda x: np.argwhere(label_values==x)[0,0])
true_labels_test=df_rad_test["IDs"].to_numpy();
label_values_test=np.unique(true_labels_test)
label_values_test.sort()

pred_hybridHMM_test = pd.read_csv('Dataset/predictions_hybrid_test.csv')
pred_hybridHMM_test['best_pep_iz'] = pred_hybridHMM_test['best_pep_iz'].apply(lambda x: np.argwhere(label_values==x)[0,0])
HMM_Hybrid_y_pred=pred_hybridHMM_test['best_pep_iz'].to_numpy();
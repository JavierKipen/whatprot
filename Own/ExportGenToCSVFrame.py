# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 17:45:32 2022

@author: JK-WORK
"""


import sys
import os
import numpy as np
#func_path="C:/Users/jkipen/Documents/GitHub/MST/Event-analysis/Data analysis/Useful_funcs_and_classes"
func_path="C:/Users/JK-WORK/Documents/whatprot/python"
if not any(os.path.normcase(sp) == os.path.normcase(func_path) for sp in sys.path):
    sys.path.append(func_path) ##Adds path to the python path to find the functions



import csv
import pandas as pd
import copy


cases=["train","test"]

for case in cases:
    ids_list=[]
    with open("Dataset/true-ids_"+case+".tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            ids_list.append(int(row[0]));
    ids_list=ids_list[1:] #Removes first label that indicates count
    
    radiometries_list=[]
    n_edman=0;n_ch=0;n_reads=0;
    with open("Dataset/radiometries_"+case+".tsv") as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for count, row in enumerate(rd):
            if count==0:
                n_edman=int(row[0])-1;
            elif count==1:
                n_ch=int(row[0]);
            elif count==2:
                n_reads=int(row[0]);
            else:
                radiometries_list.append([float(i) for i in row]);
    rad_array=np.asarray(radiometries_list)
    rads_by_char=[]
    
    for ch in range(n_ch):
        rads_by_char.append(rad_array[:,ch::n_ch]);
        
    rad_array_rearranged=copy.deepcopy(rads_by_char[0]);
    for i in range(1,n_ch):
        rad_array_rearranged=np.concatenate((rad_array_rearranged,rads_by_char[i]),axis=1)
        #rad_array_rearranged=np.insert(rad_array_rearranged, [-1], rads_by_char[i], axis=1)    
    
    t_stamps=["t_"+str(i) for i in np.arange(n_edman+1)];
    rad_stamps= ["ch_"+str(c)+"_"+t for c in range(n_ch) for t in t_stamps]
    
    final_df=pd.DataFrame(data=rad_array_rearranged, columns=rad_stamps);
    
    final_df["IDs"]=ids_list
    
    final_df.to_csv("Dataset/radiometries_"+case+".csv",index=False)


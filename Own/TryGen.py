# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:44:55 2022

@author: JK-WORK
"""

import sys
import os
import numpy as np
#func_path="C:/Users/jkipen/Documents/GitHub/MST/Event-analysis/Data analysis/Useful_funcs_and_classes"
func_path="C:/Users/JK-WORK/Documents/whatprot/python"
if not any(os.path.normcase(sp) == os.path.normcase(func_path) for sp in sys.path):
    sys.path.append(func_path) ##Adds path to the python path to find the functions


from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides

# - Replace filepaths as appropriate.
# - Replace "name-of-protease" with one of the following: "trypsin", "cyanogen bromide", or
#   "EndoPRO".
# - n is optional. If absent every protein in the .fasta file is cleaved, but if present, n
#   proteins are selected randomly without replacement from the provided file.
cleave_proteins("../examples/UP000005640_9606.fasta",
                "Dataset/peptides.csv",
                "trypsin",
                n=10)

# - Replace filepaths as appropriate.
# - label_set should be a list of strings, each string containing a set of characters
#   (typically but not always just one character). Each index of the list represents a
#   channel or color of fluorophore, and each letter represents an amino acid code.
#   For example, if we use the label_set ['DE','C','Y'], then we will label aspartic
#   acid (D) and glutamic acid (E) on channel 0, cysteine (C) on channel 1, and tyrosine
#   (Y) on channel 2.
dye_seqs_from_peptides("Dataset/peptides.csv",
                       ['DE','C','Y'],
                       "Dataset/dyeSeqs.tsv")
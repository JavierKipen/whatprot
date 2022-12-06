# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 16:44:55 2022

@author: JK-WORK
"""

from cleave_proteins import cleave_proteins
from dye_seqs_from_peptides import dye_seqs_from_peptides

# - Replace filepaths as appropriate.
# - Replace "name-of-protease" with one of the following: "trypsin", "cyanogen bromide", or
#   "EndoPRO".
# - n is optional. If absent every protein in the .fasta file is cleaved, but if present, n
#   proteins are selected randomly without replacement from the provided file.
cleave_proteins("../examples/UP000005640_9606.fasta",
                "../examples/outputfile.csv",
                "trypsin",
                n=20)

# - Replace filepaths as appropriate.
# - label_set should be a list of strings, each string containing a set of characters
#   (typically but not always just one character). Each index of the list represents a
#   channel or color of fluorophore, and each letter represents an amino acid code.
#   For example, if we use the label_set ['DE','C','Y'], then we will label aspartic
#   acid (D) and glutamic acid (E) on channel 0, cysteine (C) on channel 1, and tyrosine
#   (Y) on channel 2.
dye_seqs_from_peptides("../examples/outputfile.csv",
                       ['DE','C','Y'],
                       "../examples/dyeSeqsGen.tsv")
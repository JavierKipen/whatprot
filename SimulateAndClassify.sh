#!/bin/bash 
cd cc_code\


N_SIM_TEST=10000
N_SIM_TRAIN=300000
N_DT_PER_SEQ=1000
./bin/release/whatprot simulate rad -t 10 -g $N_SIM_TEST -P ../Own/Dataset/seqparams_atto647n_x3.json -S ../Own/Dataset/dye-seqs.tsv -R ../Own/Dataset/radiometries_test.tsv -Y ../Own/Dataset/true-ids_test.tsv
./bin/release/whatprot simulate rad -t 10 -g $N_SIM_TRAIN -P ../Own/Dataset/seqparams_atto647n_x3.json -S ../Own/Dataset/dye-seqs.tsv -R ../Own/Dataset/radiometries_train.tsv -Y ../Own/Dataset/true-ids_train.tsv
./bin/release/whatprot simulate dt -t 10 -g $N_DT_PER_SEQ -P ../Own/Dataset/seqparams_atto647n_x3.json -S ../Own/Dataset/dye-seqs.tsv -T ../Own/Dataset/dye-tracks.tsv
./bin/release/whatprot classify hybrid -k 10000 -s 0.5 -H 1000 -p 5 -P ../examples/seqparams_atto647n_x3.json -S ../Own/Dataset/dye-seqs.tsv -T ../Own/Dataset/dye-tracks.tsv -R ../Own/Dataset/radiometries_test.tsv -Y ../Own/Dataset/predictions_hybrid_test.csv

cd ..
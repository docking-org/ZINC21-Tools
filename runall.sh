#!/usr/bin/env bash

for i in x??;
do

   source /mnt/nfs/home/devtest/anaconda3/bin/activate my-rdkit-env

   python /nfs/home/jyoung/code/fine_tranche_hlogp_scripts/rdkit_hlogp_batch.py $i 10000
done

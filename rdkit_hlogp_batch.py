from __future__ import print_function
import pdb
import itertools
import sys
import shlex
import subprocess

from rdkit.Chem import (
            MolFromSmiles,
            MolToSmiles,
            AddHs,
)

from rdkit.Chem.SaltRemover import (
            SaltRemover
)

from rdkit.Chem.Descriptors import MolLogP, MolWt

def scale_logp_value(logp):
    if logp < -9.0:
        logp = -9.0
    elif logp > 9.0:
        logp = 9.0
    elif logp < 0 or logp >= 5.0:
        logp = 100*int(logp)
    elif logp >= 0.0 or logp < 1:
        logp = 10*int(10*logp)
    elif logp >= 1.0 or logp < 4:
        logp = int(100*logp)
    elif logp >=4 or logp < 5:
        logp = 10*int(10*logp)
    return logp

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python rdkit_hlogp_batch.py <smiles> <batch_size>')
        exit()
    
    BATCH_SIZE = int(sys.argv[2])
    hlogp_list = list()
    with open(sys.argv[1]) as smiles_file:
        file_lines = smiles_file.readlines()
        for line in file_lines:
            if line.strip():
                smiles, cid = str(line).strip().split()[:2]
                mol = MolFromSmiles(smiles)
                remover = SaltRemover()
                res, deleted = remover.StripMolWithDeleted(mol)
                if res is not None:
                    res.SetProp('_Name', cid)
                logp = MolLogP(res)
                num_heavy_atoms = res.GetNumHeavyAtoms()
                if num_heavy_atoms > 99:
                    num_heavy_atoms = 99
                scaled_logp = scale_logp_value(logp)
                if logp < 0.0:
                    sign = 'M'
                    #remove the minus sign so it's not printed
                    scaled_logp = scaled_logp * -1
                else:
                    sign = 'P'
                key_string = 'H{:02}{}{:03}'.format(num_heavy_atoms, sign, scaled_logp)

                #store in list up to batch size, then write out to new file
                final_string = '{0} {1} {2}\n'.format(smiles, cid, key_string)
                hlogp_list.append(final_string)

                #write the key string to the file
                if len(hlogp_list) >= BATCH_SIZE:
                    with open('{0}_hlogp'.format(sys.argv[1]), 'a+') as category_file:
                        for entry in hlogp_list:
                            category_file.write('{0}'.format(entry))
                    #clear list because we already wrote these
                    hlogp_list = list()
    #write remaining molecules to file
    for entry in hlogp_list:
        with open('{0}_hlogp'.format(sys.argv[1]), 'a+') as category_file:
                for entry in hlogp_list:
                    category_file.write('{0}'.format(entry))

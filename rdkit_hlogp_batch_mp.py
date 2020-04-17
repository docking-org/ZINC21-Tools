import multiprocessing as mp
import sys
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.SaltRemover import SaltRemover
from tqdm import tqdm

remover = SaltRemover()


def scale_logp_value(logp):
    if logp < -9.0:
        logp = -9.0
    elif logp > 9.0:
        logp = 9.0
    if logp < 0 or logp >= 5.0:
        logp = 100*int(logp)
        # logp = int(round(logp*100, -2))
    else:
        logp = 10*int(10*logp)
        # logp = int(round(logp*100, -1))
    # elif 0.0 <= logp < 1:
    #     logp = 10*int(10*logp)
    #     # logp = int(round(logp*100, -1))
    # elif 1.0 <= logp < 4:
    #     logp = int(100*logp)
    #     # logp = round(logp*100)
    # elif 4 <= logp < 5:
    #     logp = 10*int(10*logp)
    #     # logp = int(round(logp*100, -1))
    return logp


def worker(line):
    smiles, cid = line.strip().split()[:2]
    mol = MolFromSmiles(smiles)
    if mol:
        if '.' in smiles:
            mol = remover.StripMol(mol)
        logp = MolLogP(mol)
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        if num_heavy_atoms > 99:
            num_heavy_atoms = 99
        sign = 'M' if logp < 0.0 else 'P'
        return f'{smiles} {cid} H{num_heavy_atoms:02}{sign}{abs(scale_logp_value(logp)):03}\n'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python rdkit_hlogp_batch_mp.py <smiles>')
        exit()
    with open(sys.argv[1]) as f, open(f'{sys.argv[1]}_hlogp', 'w', newline='\n') as fout, mp.Pool() as pool:
        for res in tqdm(pool.imap(worker, f, chunksize=256), unit_scale=True, mininterval=1.0):
            if res:
                fout.write(res)
    print('Done')
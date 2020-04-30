import multiprocessing as mp
import sys
import os
import shutil
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
    if logp < 0.0 or logp >= 5.0:
        logp = 100*int(logp)
    else:
        logp = 10*int(10*logp)
    return logp


def worker(line):
    smiles, cid = line.decode().strip().split()[:2]
    mol = MolFromSmiles(smiles)
    if mol:
        if '.' in smiles:
            mol = remover.StripMol(mol)
        logp = MolLogP(mol)
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        if num_heavy_atoms > 99:
            num_heavy_atoms = 99
        sign = 'M' if logp < 0.0 else 'P'
        return f'{smiles} {cid}\n', len(line), True,  f'H{num_heavy_atoms:02}{sign}{abs(scale_logp_value(logp)):03}.txt'
    else:
        return f'{smiles} {cid}\n', len(line), False


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python rdkit_hlogp_batch_mp_2.py <smiles>')
        exit()
    folder = os.path.dirname(sys.argv[1])
    name, ext = os.path.splitext(os.path.basename(sys.argv[1]))
    file_failed = os.path.join(folder, f'{name}_failed{ext}')
    tranches_folder = f'./{name}_tranches/'
    if os.path.exists(tranches_folder):
        shutil.rmtree(tranches_folder)
    os.makedirs(tranches_folder)
    cache = dict()
    cache_size_in_lines = 1024
    with open(sys.argv[1], 'rb') as f, open(file_failed, 'w', newline='\n') as f_f, mp.Pool() as pool, \
            tqdm(total=os.path.getsize(sys.argv[1]), unit_scale=True, unit_divisor=1024,
                 unit='B', mininterval=1.0) as _tqdm:
        for res in pool.imap(worker, f, chunksize=256):
            _tqdm.update(res[1])
            if res[2]:
                cache.setdefault(res[3], list()).append(res[0])
                if len(cache[res[3]]) >= cache_size_in_lines:
                    with open(os.path.join(tranches_folder, res[3]), 'a', newline='\n') as fout:
                        fout.writelines(cache[res[3]])
                    cache[res[3]].clear()
            else:
                f_f.write(res[0])
    for file_name, cached_lines in cache.items():
        if cached_lines:
            with open(os.path.join(tranches_folder, file_name), 'a', newline='\n') as fout:
                fout.writelines(cached_lines)
    print('Done')

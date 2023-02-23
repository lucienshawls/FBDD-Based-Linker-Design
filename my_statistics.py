import csv
import re
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


TOTAL = 3270

def is_valid(smiles):
    
    mol = Chem.MolFromSmiles(smiles)

    mw = Descriptors.MolWt(mol)
    nc = len(mol.GetSubstructMatches(Chem.MolFromSmiles('C')))
    nar = Chem.Lipinski.NumAromaticRings(mol)
    hbd = Chem.Lipinski.NumHDonors(mol)
    hba = Chem.Lipinski.NumHAcceptors(mol)

    tpsa = Descriptors.TPSA(mol)
    # print(mw,nc,nar,hbd,hba,tpsa)
    if (10<=mw<=600) and (1<=nc<=25) and (0<=nar<=2) and (0<=hbd<=7) and (0<=hba<=15) and (0<=tpsa<=180):
        return True
    else:
        return False

def count_fragments(verification='skip',dataset_base_dir='./data/macfrag/protac/'):
    res = []
    tot = TOTAL
    for k in range(1,tot+1):
        detail = {'protac_id': k}
        count = 0
        for i in range(3,10+1):
            for j in range(5,8+1):
                para = '%d-%d-F-1' %(i,j)
                frag_path = '%s%s/protac_%d_%s_fragments.smi' %(dataset_base_dir,para,k,para)
                with open(frag_path, 'r', encoding='utf-8') as f:
                    txt = f.read().strip()
                    if verification == 'skip':
                        cnt = len(txt.split('\n'))
                    elif verification == 'linker_validation':
                        print('Counting valid linkers of protac: %d (for parameter %s)' %(k,para),end='')
                        cnt = 0
                        for linker in txt.split('\n'):
                            linker_refined = re.sub('\[[0-9]*\*\]','[*]',linker)
                            if is_valid(linker_refined):
                                cnt += 1
                        print(', total: %d' %(cnt))
                    count += cnt
                    detail[para] = cnt
        detail['TOTAL'] = count
        res.append(detail)

    detail = {'protac_id': 'TOTAL'}
    count = 0
    for i in range(3,10+1):
        for j in range(5,8+1):
            para = '%d-%d-F-1' %(i,j)
            cnt = 0
            for k in range(1,tot+1):
                cnt += res[k-1][para]
            count += cnt
            detail[para] = cnt
    detail['TOTAL'] = count
    res.append(detail)

    return res

def write_file(filename, mydict):
    with open(filename,'w',encoding='utf-8',newline='') as f:
        my_header = list(mydict[0].keys())
        writer = csv.DictWriter(f,my_header)
        writer.writeheader()
        writer.writerows(mydict)

def main():
    # count = count_fragments()
    # write_file('./data/results/total.csv',count)
    count = count_fragments('linker_validation')
    write_file('./data/results/total_valid_linker.csv',count)
    # a = is_valid('*c1ccccc1O*')
    # print(a)
    pass


if __name__ == '__main__':
    main()

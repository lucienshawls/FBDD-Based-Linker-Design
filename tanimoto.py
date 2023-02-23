import csv
import random
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import re

TOTAL = 3270
SIMILARITY_THRESHOLD = 0.7

def similarity(reference,test):
    # print('{%s},{%s}' %(reference,test))
    # res = random.random()
    mol_r = Chem.MolFromSmiles(reference)
    mol_t = Chem.MolFromSmiles(test)
    fps_r = Chem.RDKFingerprint(mol_r)
    fps_t = Chem.RDKFingerprint(mol_t)

    # fps_r = AllChem.GetMorganFingerprintAsBitVect(mol_r,2)
    # fps_t = AllChem.GetMorganFingerprintAsBitVect(mol_t,2)
    # res = DataStructs.FingerprintSimilarity(fps_r,fps_t)
    res = DataStructs.DiceSimilarity(fps_r,fps_t)
    return res

def recall_rate(para,
                reference_path_protac='./data/protac.csv',
                reference_path_linker='./data/linker.csv',
                test_database_dir='./data/macfrag/protac/'):
    tot = TOTAL
    positive = 0
    invalid = 0
    negative = 0
    for i in range(1,tot+1):
        print('Calculating recall rate: %d (for parameter %s)' %(i,para),end='')
        with open(reference_path_protac, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            linker_id = list(reader)[i-1]['linker_id'] # id starts from 0, not 1
            # reader = csv.reader(f)
            # reference = list(reader)[i][6]
        if linker_id == '':
            print(', NO_LINKER')
            invalid += 1
            continue
        with open(reference_path_linker,'r',encoding='utf-8') as f:
            reader = csv.DictReader(f)
            reference = list(reader)[int(linker_id)-1]['Smiles_R'] # id starts from 0, not 1
            reference_refined = re.sub('\[R[0-9]*\]','[*]',reference.strip())
        test_database_path = '%s%s/protac_%d_%s_fragments.smi' %(test_database_dir,para,i,para)
        test_database = []
        with open(test_database_path, 'r', encoding='utf-8') as f:
            test_database = list(f.readlines())
        is_match = False
        for test in test_database:
            test_refined = re.sub('\[[0-9]*\*\]','[*]',test.strip())
            if similarity(reference_refined,test_refined) >= SIMILARITY_THRESHOLD:
                is_match = True
                break
        if is_match:
            print(', match_found: yes, {%s} vs {%s}' %(reference_refined, test_refined))
            positive += 1
        else:
            print(', match_found: no')
            negative += 1
    print('positive: %d\n' %(positive))
    rate = positive/tot
    return rate

def para_assess():
    res = []
    
    with open('./data/results/recall_rate.csv', 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['parameter','recall_rate'])

    for i in range(3,10+1):
        for j in range(5,8+1):
            para = '%d-%d-F-1' %(i,j)
            record = {
                'parameter': para,
                'recall_rate': recall_rate(para)
            }
            with open('./data/results/recall_rate.csv', 'a', encoding='utf-8', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([record['parameter'],record['recall_rate']])
            res.append(record)
    return res

def main():
    recall_rate_set = para_assess()
    pass

if __name__ == '__main__':
    main()

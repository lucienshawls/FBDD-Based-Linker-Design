import csv


TOTAL = 3270
SIMILARITY_THRESHOLD = 0.9

def similarity(reference,test):
    res = 1
    return res

def recall_rate(para,reference_path='./data/protac.csv',test_database_dir='./data/macfrag/protac/'):
    tot = TOTAL
    positive = 0
    for i in range(1,tot+1):
        with open(reference_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            reference = list(reader)[i-1]['linker_smiles'] # id starts from 0, not 1

        test_database_path = '%s%s/protac_%d_%s_fragments.smi' %(test_database_dir,para,i,para)
        test_database = []
        with open(test_database_path, 'r', encoding='utf-8') as f:
            for row in f.readlines():
                test_database.append(str(row))
        is_match = False
        for test in test_database:
            if similarity(reference,test) >= SIMILARITY_THRESHOLD:
                is_match = True
                break
        if is_match:
            positive += 1

    rate = positive/tot
    return rate

def para_assess():
    res = []
    for i in range(3,10+1):
        for j in range(5,8+1):
            para = '%d-%d-F-1' %(i,j)
            res.append({
                'parameter': para,
                'recall_rate': recall_rate(para)
            })
    return res
            
def main():
    recall_rate_set = para_assess()
    with open('./data/results/recall_rate.csv', 'w', encoding='utf-8') as f:
        writer = csv.DictWriter(f,['parameter','recall_rate'])
        writer.writeheader()
        writer.writerows(recall_rate_set)

if __name__ == '__main__':
    main()

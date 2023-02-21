import requests
import csv

COMPOUND_LIST_PAGE = 'http://cadd.zju.edu.cn/protacdb/browse/compound'
TOTAL = 3270

def compound(id):
    info = {
        'id': id,
        'target': '',
        'protac': {
            'id': '',
            'smiles': ''
        },
        'e3_ligand': {
            'id': '',
            'smiles': ''
        },
        'linker': {
            'id': '',
            'smiles': ''
        },
        'warhead': {
            'id': '',
            'smiles': ''
        },
    }
    compound_page = 'http://cadd.zju.edu.cn/protacdb/compound/dataset=protac&id=%d' %(id)

    cookies = {'_xsrf': '2|77ac8144|82b88bc11af2f09bfdf53c753ca8b906|1676875379',}
    headers = {
        'Accept': '*/*',
        'Accept-Language': 'en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7,en-GB;q=0.6',
        'Connection': 'keep-alive',
        'Content-Type': 'multipart/form-data; boundary=----WebKitFormBoundaryRJRMQVyUyjFMx6O0',
        'DNT': '1',
        'Origin': 'http://cadd.zju.edu.cn',
        'Referer': compound_page,
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36 Edg/110.0.1587.50',
        'X-Requested-With': 'XMLHttpRequest',
    }
    data = '------WebKitFormBoundaryRJRMQVyUyjFMx6O0\r\nContent-Disposition: form-data; name="_xsrf"\r\n\r\n2|77ac8144|82b88bc11af2f09bfdf53c753ca8b906|1676875379\r\n------WebKitFormBoundaryRJRMQVyUyjFMx6O0--\r\n'

    response = requests.post(
        compound_page,
        headers=headers,
        verify=False,
        cookies=cookies,
        data=data
    )
    resp = response.json()

    info['protac']['id'] = resp[0]['id_in_database']
    info['protac']['smiles'] = resp[0]['smiles_canonical']
    info['e3_ligand']['id'] = resp[0]['id_e3_ligand']
    info['e3_ligand']['smiles'] = resp[0]['e3_ligand_canonical']
    info['linker']['id'] = resp[0]['id_linker']
    info['linker']['smiles'] = resp[0]['linker_canonical']
    info['warhead']['id'] = resp[0]['id_warhead']
    info['warhead']['smiles'] = resp[0]['warhead_canonical']

    info['target'] = resp[-1]['target']
    return info

def compound_list():
    tot = TOTAL
    dataset = []

    with open('./data/protac.csv','w',encoding='utf-8-sig',newline='') as f:
        writer = csv.writer(f)

        cur_row = [
            'id',
            'protac_id',
            'protac_smiles',
            'e3_ligand_id',
            'e3_ligand_smiles',
            'linker_id',
            'linker_smiles',
            'warhead_id',
            'warhead_smiles',
            'target',
        ]
        writer.writerow(cur_row)

    for i in range (1,tot+1):
        print('Processing: %d/%d' %(i,tot))
        data = compound(i)
        cur_row = [
            data['id'],
            data['protac']['id'],
            data['protac']['smiles'],
            data['e3_ligand']['id'],
            data['e3_ligand']['smiles'],
            data['linker']['id'],
            data['linker']['smiles'],
            data['warhead']['id'],
            data['warhead']['smiles'],
            data['target'],
        ]
        with open('./data/protac.csv','a',encoding='utf-8-sig',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(cur_row)
        dataset.append(data)
    
    return dataset

def main():
    dataset = compound_list()

if __name__ == '__main__':
    main()

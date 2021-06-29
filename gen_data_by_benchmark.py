# -*- coding: utf-8 -*-
# @Time       : 2021/06/23 10:34:55
# @Author     : Yuang Tong <yuangtong1999@gmail.com>
# @Project    : CPI-GLX4.1.0
# @Description:  根据benchmark建立数据集

'''
根据benchmark的蛋白质，找到600w中对应的蛋白，构建另一个更小的数据集
'''

import pandas as pd
import os
import numpy as np
from tqdm import tqdm

raw_data = pd.read_csv('/gxr/tongyuang/data/CPI/GLX4.1.0/filtered_data/Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7_rmv_dup.csv')

benchmark_Data = pd.read_csv('/gxr/siqi/8000w_data_folder/seen_protein_non_structure_benchmark_v2.csv')

bm_sqs_list = (list(set(list(benchmark_Data['PROTEIN_SEQUENCE'].values))))


COL_NAMES=['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']
PRT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/'
DATA_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx4.1.0_clf_Pubchem_Uniprot_Classification_Data_79980610/'
SQS_IDX_DICT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/sequence_idx_dict.npy'
sqs_idx_dict = np.load(SQS_IDX_DICT_ROOT,allow_pickle=True).item()

def index_to_name(subdata_index):
    '''
    input: int
        0->00
        19->19
    output: str
    '''
    if subdata_index<10:
        return '0'+str(int(subdata_index))
    else:
        return str(int(subdata_index))

def sqs_subdata(sqs_in='',
                sqs_idx_dict=sqs_idx_dict,
                prt_root= PRT_ROOT,
                data_root = DATA_ROOT,
                col_names = COL_NAMES,
                suffix='.npy'):
    '''
    get subdf by sequence

    输入：蛋白质sequence
    输出：subdf
    '''

    subdf = pd.DataFrame(columns=col_names)

    prt_idx = sqs_idx_dict[sqs_in]
    filename = str(prt_idx)+suffix
    cur_prt_statistics = np.load(prt_root+filename,allow_pickle=True)
    
    subdata_index_list = list(set(cur_prt_statistics[:,0]))
    
    output_smiles_list = []
    
    for subdata_index in subdata_index_list:
        
        subdata_index = int(subdata_index)
        
        rows = cur_prt_statistics[cur_prt_statistics[:,0]==subdata_index][:,1]
        subdata_name = index_to_name(subdata_index)
        

        if subdata_name=='00':
            cur_df = pd.read_csv(data_root+'part_{}.csv'.format(subdata_name))
        else:
            cur_df = pd.read_csv(data_root+'part_{}.csv'.format(subdata_name),names=col_names)
            
        sub_cur_df = cur_df.iloc[rows]
        

        subdf = subdf.append(sub_cur_df)
        #smiles_list = list(sub_cur_df['COMPOUND_SMILES'].values)
            
        #output_smiles_list += smiles_list

    return subdf

output_df = pd.DataFrame(columns=COL_NAMES)

for sqs in tqdm(bm_sqs_list):
    # subdata
    try:
        subdf = sqs_subdata(sqs_in=sqs)
        # remove duplicate
        ref_smiles = list(benchmark_Data[benchmark_Data['PROTEIN_SEQUENCE']==sqs]['COMPOUND_SMILES'].values)
        drop_df = subdf[subdf['COMPOUND_SMILES'].isin(ref_smiles)]
        subdf = subdf.drop(index=list(drop_df.index))

        # append
        output_df = output_df.append(subdf,ignore_index=True)
        #break
    
    except:

        print('{} is not in the dict'.format(sqs))
        continue

print('saving...')
output_df.to_csv('/gxr/tongyuang/data/CPI/GLX4.1.0/filtered_data/Glx_Pubchem_Uniprot_Classification_Data_raw_from_bm.csv',index=0)
print('Done.')


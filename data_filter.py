# -*- coding: utf-8 -*-
# @Time       : 2021/06/16 17:16:33
# @Author     : Tong Yuang <yuangtong1999@gmail.com>
# @Project    : CPI Data Analysis
# @Description: 过滤8kw数据

import pandas as pd
import numpy as np
from sklearn.utils import shuffle
import os

from threading import Thread # 并行处理


COL_NAMES=['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']
PRT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/'
DATA_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx4.1.0_clf_Pubchem_Uniprot_Classification_Data_79980610/'
SQS_IDX_DICT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/sequence_idx_dict.npy'
OUTPUT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/filted/Glx_Pubchem_Uniprot_Classification_Data.csv'
RESULT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/filted/Glx_Pubchem_Uniprot_Classification_Data'

sqs_idx_dict = np.load(SQS_IDX_DICT_ROOT,allow_pickle=True).item()
statistics_file_name = './sequence_statistics.csv'
statistics_df = pd.read_csv(statistics_file_name)

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



from tqdm import tqdm

def filt_data(index_range,thres=0.18, res_root = RESULT_ROOT,columns=COL_NAMES,ref_df = statistics_df):

    '''
    index_range:[(0,100),(1000,1100),...]

    '''
    for idx_pair in index_range:

        
        start = idx_pair[0]
        end = idx_pair[1]
        print('filting range:{}-{}'.format(start,end))

        cur_root = res_root+'_{}.csv'.format(start)

        if not os.path.exists(cur_root):
            output_df = pd.DataFrame(columns=columns)
            output_df.to_csv(cur_root)

        output_df = pd.read_csv(cur_root)

        for i in tqdm(range(start,end)):

            if i%20 == 19:
                output_df.to_csv(cur_root)
    
            row = ref_df.iloc[i]
            sqs = row['sequence']
            pos_pct = row['pos_percent']

            subdf = sqs_subdata(sqs_in=sqs)

            if pos_pct<thres:
                new_subdf = pd.DataFrame(columns=columns)

                pos_subdf = subdf[subdf['CLF_LABEL']==1]
                neg_subdf = subdf[subdf['CLF_LABEL']==0]

                target_num = (int)(len(pos_subdf)/thres)+1

                del_num = len(neg_subdf)-target_num

                neg_subdf = shuffle(neg_subdf)

                neg_subdf.drop(neg_subdf.head(del_num).index,inplace=True)

                new_subdf = new_subdf.append(pos_subdf)
                new_subdf = new_subdf.append(neg_subdf)
        
                new_subdf = shuffle(new_subdf)

                output_df = output_df.append(new_subdf)
            else:
                output_df = output_df.append(subdf)
    
        output_df.to_csv(cur_root)

    print('Done')
    return None

def main():


    num_threads = 10
    raw_index = np.arange(len(statistics_df))
    interval = 100
    start = 0
    end = start+interval-1

    index_range_list = []

    for i in range(num_threads):
        index_range_list.append([])

    while(True):

        cur_row = (int)(start/interval)%num_threads

        index_range_list[cur_row].append((start,end))

        start = end+1
        end = start+interval-1
        if end>=len(statistics_df)-1:
            end = len(statistics_df)-1
            cur_row = (int)(start/interval)%num_threads
            index_range_list[cur_row].append((start,end))
            break


    threads = []
    for i in range(num_threads):
        thread = Thread(target=filt_data,args=(index_range_list[i],0.18, RESULT_ROOT,COL_NAMES,statistics_df))

        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()


    print('Done')

if __name__ == '__main__':
    main()
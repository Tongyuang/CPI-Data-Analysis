# created by Yuang Tong at 2021.6.3

import pandas as pd
import numpy as np
from RdKit_tools import calc_fingerprints,calc_similarity_matrix
from tqdm import tqdm
import time

import random
import os

random.seed(2021)

COL_NAMES=['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']
PRT_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/'
DATA_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx4.1.0_clf_Pubchem_Uniprot_Classification_Data_79980610/'

output_file_name = './sequence_statistics.csv'
output_df = pd.read_csv(output_file_name)


def gen_similarity_matrix(smlies_list):
    '''
    input: smiles_list
    output: similarity_matrix
    '''
    return calc_similarity_matrix(calc_fingerprints(smlies_list))


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

def get_smiles_list(prt_idx=0,
                    prt_root= PRT_ROOT,
                    data_root = DATA_ROOT,
                    col_names = COL_NAMES,
                   suffix='.npy'):
    '''
    get smiles list by protein index
    '''
    filename = str(prt_idx)+suffix
    cur_prt_statistics = np.load(prt_root+filename,allow_pickle=True)
    
    subdata_index_list = list(set(cur_prt_statistics[:,0]))
    
    output_smiles_list = []
    
    for subdata_index in tqdm(subdata_index_list):
        
        subdata_index = int(subdata_index)
        
        rows = cur_prt_statistics[cur_prt_statistics[:,0]==subdata_index][:,1]
        subdata_name = index_to_name(subdata_index)
        

        if subdata_name=='00':
            cur_df = pd.read_csv(data_root+'part_{}.csv'.format(subdata_name))
        else:
            cur_df = pd.read_csv(data_root+'part_{}.csv'.format(subdata_name),names=col_names)
            
        sub_cur_df = cur_df.iloc[rows]
            
        smiles_list = list(sub_cur_df['COMPOUND_SMILES'].values)
            
        output_smiles_list += smiles_list

    return list(set(output_smiles_list))
            

def cal_silimarity_metrics(smiles_list):
    
    if(len(smiles_list)==1):
        print('only one compound smile record.')
        return {'mean_sim':1,
        'max_sim':1,
        'var_sim':0
    }
    if(len(smiles_list)>1000):
        update_smile_list = random.sample(smiles_list, 1000)
        print('input smile list is too long (>{}), random choose {} for evaluation'.format(1000,1000))
    else:
        update_smile_list = smiles_list
    fplist = calc_fingerprints(update_smile_list)

    sim_matrix = calc_similarity_matrix(fplist)
    #print(len(sim_matrix))
    sim_triup = np.triu(sim_matrix, 1).flatten()
    #print(sim_triup[0:10])
    sim_list = []
    for value in sim_triup:
        if value>0:
            sim_list.append(value)
    
    if len(sim_list)==0:
        print('similarity=0')
        return {'mean_sim':0,
                'max_sim':0,
                'var_sim':0
                }

    return {'mean_sim':np.mean(sim_list),
            'max_sim':np.max(sim_list),
            'var_sim':np.var(sim_list)
    }

prt_index_list = list(output_df['idx'].values)

tic = time.time()
timeArray = time.localtime(tic)




start_point = 945


for i, idx in enumerate(prt_index_list[start_point:],start=start_point):

    smiles_list = get_smiles_list(prt_idx=idx)
    
    sim_dict = cal_silimarity_metrics(smiles_list)
    
    output_df.iloc[i,6] = sim_dict['mean_sim']
    output_df.iloc[i,7] = sim_dict['max_sim']
    output_df.iloc[i,8] = sim_dict['var_sim']
    
    if i%10 == 9:
        print('finished {:d} / {:d}'.format(i+1,len(prt_index_list)))
        toc = time.time()
        print('Time cost:{:.2f} s'.format(toc-tic))
        output_df.to_csv(output_file_name,index=False)
        tic = time.time()

print('Done.')
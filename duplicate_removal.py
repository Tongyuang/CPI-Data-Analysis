# -*- coding: utf-8 -*-
# @Time       : 2021/06/17 16:13:07
# @Author     : Yuang Tong <yuangtong1999@gmail.com>
# @Project    : data_analysis
# @Description:  数据集去重

import pandas as pd
import os
from tqdm import tqdm


COL_NAMES=['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']
SUBDF_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/filted/'
output_name = 'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7.csv'
OUTPUT_ROOT= '/gxr/tongyuang/data/CPI/GLX4.1.0/filtered_data/'

'''
output_df = pd.DataFrame(columns=COL_NAMES)
output_name = 'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7_rmv_dup.csv'
output_df.to_csv(OUTPUT_ROOT+output_name,index=0)
'''
output_df=pd.read_csv(OUTPUT_ROOT+'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7_rmv_dup.csv')

benchmark_df = pd.read_csv('/gxr/siqi/8000w_data_folder/seen_protein_non_structure_benchmark_v2.csv')

bm_sqs_list = list(set(list(benchmark_df['PROTEIN_SEQUENCE'].values)))


for subdf_name in tqdm(os.listdir(SUBDF_ROOT)):
    subdf = pd.read_csv(SUBDF_ROOT+subdf_name)
    subdf.drop(['Unnamed: 0','Unnamed: 0.1'],axis=1,inplace=True)

    temp = subdf[subdf['PROTEIN_SEQUENCE'].isin(bm_sqs_list)]
    temp_sqs_list = list(set(list(temp['PROTEIN_SEQUENCE'].values)))



    drop_df_all = pd.DataFrame(columns=temp.columns)

    for sqs in temp_sqs_list:

        target = temp[temp['PROTEIN_SEQUENCE']==sqs]

        ref = benchmark_df[benchmark_df['PROTEIN_SEQUENCE']==sqs]

        drop_df = target[target['COMPOUND_SMILES'].isin(list(ref['COMPOUND_SMILES'].values))]

        drop_df_all = drop_df_all.append(drop_df)
    #print(subdf.columns)
    subdf = subdf.drop(index=list(drop_df_all.index))

    output_df = output_df.append(subdf,ignore_index=True)


print('saving...')
output_name = 'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7_rmv_dup.csv'
output_df.to_csv(OUTPUT_ROOT+output_name,index=0)
print('Done.')
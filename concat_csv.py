# -*- coding: utf-8 -*-
# @Time       : 2021/06/17 10:25:05
# @Author     : Yuang Tong <yuangtong1999@gmail.com>
# @Project    : data_analysis
# @Description:  把小数据集拼接成大数据集


import pandas as pd
import os
from tqdm import tqdm
COL_NAMES=['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']
SUBDF_ROOT = '/gxr/tongyuang/data/CPI/GLX4.1.0/filted/'
output_name = 'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7.csv'
OUTPUT_ROOT= '/gxr/tongyuang/data/CPI/GLX4.1.0/filtered_data/'

output_df=pd.read_csv(OUTPUT_ROOT+'Glx_Pubchem_Uniprot_Classification_Data_DROP_1_7.csv')

for subdf_name in tqdm(os.listdir(SUBDF_ROOT)):
    subdf = pd.read_csv(SUBDF_ROOT+subdf_name)
    subdf.drop(['Unnamed: 0','Unnamed: 0.1'],axis=1,inplace=True)
    #print(subdf.columns)
    output_df = output_df.append(subdf,ignore_index=True)

print('saving...')
output_df.to_csv(OUTPUT_ROOT+output_name,index=0)
print('Done.')
# created by Yuang Tong at 2021.5.26

import os
import time
import numpy as np
import cudf

protein_names = np.load('protein_array.npy',allow_pickle=True)
#data_path = '/gxr/yize/data/CPI/GLX4.1.0/CLF/Glx_Pubchem_Uniprot_Classification_Data_79980610.csv'
#data_path = '/gxr/yize/data/CPI/GLX4.1.0/ALL/Glx4.1.0_ALL_IC50.csv'
os.environ['CUDA_VISIBLE_DEVICES']='4'
col_names = ['UID','UNIPROT_ID','COMPOUND_SMILES','PROTEIN_SEQUENCE','CLF_LABEL','Data_Source']

data_root = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx4.1.0_clf_Pubchem_Uniprot_Classification_Data_79980610'
save_root = '/gxr/tongyuang/data/CPI/GLX4.1.0/sub_data/Glx_Pubchem_Uniprot_Classification_Data_79980610_proteins/'

extension = '.csv'
file_name_list = []
for filename in os.listdir(data_root):
    if filename.endswith(extension):
        file_name_list.append(data_root+'/'+filename)
file_name_list.sort()

prt_dict = {}
for name in protein_names:
    prt_dict[name] = []

sequence_idx_dict = {}
for i in range(len(protein_names)):
    sequence_idx_dict[protein_names[i]] = i
np.save(save_root+'sequence_idx_dict.npy',sequence_idx_dict)

print('sequence_idx_dict saved at {}'.format(save_root))

for i in range(len(file_name_list)):
    if(i==0):
        df = cudf.read_csv(file_name_list[i])
    else:
        df = cudf.read_csv(file_name_list[i],names=col_names)
    
    print('subdata {} loaded, subdata length={}'.format(file_name_list[i],len(df)))
    
    print('iterating subdata {}/{}'.format(i,len(file_name_list)-1))
    
    subdata_idx = i
    
    tic = time.time()
    timeArray = time.localtime(tic)
    print('start iterating:'+time.strftime("%Y-%m-%d %H:%M:%S", timeArray))
    for idx,row in df.to_pandas().iterrows():
        sequence = row['PROTEIN_SEQUENCE']
        subrow_idx = idx
        row_idx = subdata_idx*(1000000)+subrow_idx
        
        #smiles = row['COMPOUND_SMILES']
        label = row['CLF_LABEL']
        
        prt_dict[sequence].append([subdata_idx,subrow_idx,row_idx,label])
    toc = time.time()
    timeArray = time.localtime(toc)
    print('end iterating:'+time.strftime("%Y-%m-%d %H:%M:%S", timeArray)+'. Time cost:{:.2f} s'.format(toc-tic))
    #break
print('End iterating.')

# save
print('saving to {}'.format(save_root))
for key in prt_dict.keys():
    # name_in_list = prt_dict[key]
    
    name = sequence_idx_dict[key]
    print('writing protein {}'.format(name))
    np.save(save_root+str(name)+'.npy',prt_dict[key])
    
    #break
print('Done.')
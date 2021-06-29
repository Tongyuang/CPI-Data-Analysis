# -*- coding: utf-8 -*-
# @Time       : 2021/03/31 09:50:33
# @Author     : Jiang Yize <yize.jiang@galixir.com>
# @Project    : Galixir-Comet-Ranking
# @Description: RDKit相关工具

from abc import abstractmethod
import logging
from typing import List, Tuple, Union
from multiprocessing import Pool

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm


pool = Pool(10)

def calc_fingerprints_values(fp_in):
    '''
    计算fingerprint 的值
    '''
    return [float(i) for i in fp_in.ToBitString()]

def calc_fingerprints(smiles: Union[List[str], pd.DataFrame]):

    """计算分子的Fingerprint特征列表

    Args:
        smiles (List[str]): 分子的SMILES列表

    Returns:
        List: 分子的Fingerprint列表"""
    
    if isinstance(smiles, pd.DataFrame):
        smiles = smiles['COMPOUND_SMILES']
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in mols]
    return fps


def calc_fingerprints_multi_threads(smiles: Union[List[str], pd.DataFrame]):
    """
    多进程计算分子的Fingerprint特征列表

    Args:
        smiles (List[str]): 分子的SMILES列表

    Returns:
        List: 分子的Fingerprint列表
    """
    if isinstance(smiles, pd.DataFrame):
        smiles = smiles['COMPOUND_SMILES']
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    
    def get_fp(x):
        return AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024)

    fps = []
    for each in tqdm(pool.imap(get_fp, mols), desc='calculating fps'):
        fps.append(each)
    
    return fps


def calc_dist_matrix(input_data: pd.DataFrame):
    """计算Tanimoto距离矩阵

    Args:
        input_data (pd.DataFrame): 输入数据

    Returns:
        List: Tanimoto距离矩阵
    """

    logging.info("\tStart Calculate Distance Matrix ... ")
    fps = calc_fingerprints(input_data)
    dist_mat = [DataStructs.BulkTanimotoSimilarity(fps[i], fps, returnDistance=True) for i in range(len(fps))]
    logging.info("\tFinish Calculate Distance Matrix.")
    return dist_mat


def calc_distance_matrix(fps: List[int]):
    """计算Tanimoto距离矩阵（三角矩阵）

    Args:
        fps (List[int]): 分子的Fingerprint列表

    Returns:
        List[float]: 距离矩阵（三角矩阵）
    """
    dist_mat = []
    for i in tqdm(range(1, len(fps)), desc="Calculate Distance Matrix"):
        simi = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dist_mat.extend([1 - x for x in simi])
    return dist_mat


def calc_similarity_matrix(fps: List[int]):
    """计算Tanimoto距离相似度矩阵

    Args:
        fps (List[int]): 分子的Fingerprint列表

    Returns:
        List[float]: 相似度矩阵
    """
    simi_mat = [DataStructs.BulkTanimotoSimilarity(fps[i], fps, returnDistance=False) for i in range(len(fps))]
    return simi_mat


def calc_compound_similarity(smi1: str, smi2: str):
    """计算两个分子的相似度

    Args:
        smi1 (str): 分子的SMILES
        smi2 (str): 分子的SMILES

    Returns:
        float: 相似度
    """
    fp1 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi1), 2, 1024)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi2), 2, 1024)
    simi = DataStructs.TanimotoSimilarity(fp1, fp2)
    return simi


def calc_reference_similarity(ref_smiles_list, curr_smiles):

    """计算某一分子与多个参考分子的相似度

    Args:
        ref_smiles_list ([type]): 参考分子列表
        curr_smiles ([type]): 当前分子

    Returns:
        Tuple[str, float]: 最相似的分子和与其的相似度
    """
    ref_fps = []
    for s in ref_smiles_list:
        ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 1024))
    
    curr_fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(curr_smiles), 2, 1024)
    simi = DataStructs.BulkTanimotoSimilarity(curr_fp, ref_fps)
    index = np.argmax(simi)
    return ref_smiles_list[index], np.max(simi)


def get_similarity_matrix(ref_data, commercial_data, save_path=None):
    
    """

    Args:
        ref_data: 参考分子
        commercial_data: 商业库分子
        save_path: 相似度矩阵文件存储的位置

    Returns:
        np.ndarray: 相似度矩阵，维度为(n_ref, n_commercial)
    """
    ref_smiles = ref_data["COMPOUND_SMILES"].tolist()

    ref_fps = []
    for s in ref_smiles:
        ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 1024))

    comm_fps = []
    for s in tqdm(commercial_data["COMPOUND_SMILES"].tolist(), desc='Calculating commercial data fp...'):
        comm_fps.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 1024))
    #comm_fps = calc_fingerprints_multi_threads(commercial_data["COMPOUND_SMILES"].tolist())
    simi_matrix = np.zeros((len(ref_fps), len(comm_fps)))
    for i in tqdm(range(len(ref_fps)), desc='Calculating similarity matrix...'):
        simi_matrix[i] = DataStructs.BulkTanimotoSimilarity(ref_fps[i], comm_fps)

    if save_path is not None:
        np.save(save_path, simi_matrix)
        logging.info(f"Save similarity matrix to {save_path}")
    return simi_matrix


def filter_cluster_data(
    cluster_data: pd.DataFrame,
    filter_topn: int,
    rank_column: str,
    use_murcko: bool = True,
    cluster_simi_threshold: float = 0.7
):
    """对聚类后的数据进行过滤

    Args:
        cluster_data (pd.DataFrame): 包含ClusterID的数据
        filter_topn (int): 每个类簇选取的TopN
        rank_column (str): 排序依赖的列名
        use_murcko (bool, optional): 是否使用骨架ID过滤. Defaults to True.
        cluster_simi_threshold:
    Returns:
        pd.DataFrame: 聚类过滤后的数据
    """


    n_cluster = len(set(cluster_data["ClusterID"]))

    cluster_filter_data_list = []
    for i in tqdm(range(n_cluster)):
        temp_data = cluster_data[cluster_data["ClusterID"] == i].copy()
        temp_data = temp_data.sort_values(by=rank_column, ascending=False)

        select_idxs = [0]
        for j in range(1, len(temp_data)):
            current_smiles = temp_data.iloc[j, :]["COMPOUND_SMILES"]
            for select_idx in select_idxs:
                select_smiles = temp_data.iloc[select_idx, :]["COMPOUND_SMILES"]
                if calc_compound_similarity(current_smiles, select_smiles) >= cluster_simi_threshold:
                    break
            else:
                select_idxs.append(j)
            if len(select_idxs) >= filter_topn:
                break
        cluster_filter_data_list.append(temp_data.iloc[select_idxs, :])

    cluster_filter_data = pd.concat(cluster_filter_data_list)
    cluster_filter_data = cluster_filter_data.sort_values(by=rank_column, ascending=False)
    logging.info(f"\tCluster TopN Number: {len(cluster_filter_data)}")
    if use_murcko:
        logging.info("\tUse Murcko Filter!")
        cluster_filter_data = cluster_filter_data.drop_duplicates(subset="MurckoID", keep="first")

    return cluster_filter_data


if __name__ == '__main__':
    
    
    test_smiles = [
        'Brc1ccc(cc1)C1CC(=NN1C(=O)C)c1ccccc1',
        'Brc1ccc(cc1)C1CC(=NN1c1ccccc1)c1ccccc1',
        'Brc1ccc(cc1)C1NN=C(C1)c1ccccc1'
        ]
    
    #fplist = calc_fingerprints(test_smiles)

    #sim_matrix = calc_similarity_matrix(fplist)

    #print(sim_matrix)

    test_fps = calc_fingerprints(test_smiles)

    for fp in test_fps:
        print(len(calc_fingerprints_values(fp)))


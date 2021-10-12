#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description:
# usage: python3 batch_GCN [specificed_atom_id_ASE_num] [CN_max]
import os
from ase.io import read, write
import pandas as pd
import numpy as np
import sys
from alive_progress import alive_bar
import time


dict_Pyykko_radii = {'O': 0.63, 'Sc': 1.49, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19,
                     'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1, 'Cu': 1.12, 'Zn': 1.18, 'Y': 1.63,
                     'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
                     'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36, 'Ce': 1.63, 'Hf': 1.52, 'Ta': 1.46,
                     'W': 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24,
                     'Hg': 1.33 }

def CN( atom_id ):
    CN_SUM = 0
    atom_id = int(atom_id)
    Ra = dict_Pyykko_radii[symbols[atom_id]]
    for i in range(len(distances_atom2others)):
            if i != int(atom_id):
                Rb = dict_Pyykko_radii[symbols[i]]
                Rab = distances_atom2others[i]
                # Rab = slab.get_distance(0, 2, mic=False, vector=False)
                CN_sub = 1 / (1 + np.exp(-16 * (4/3 * (Ra + Rb)/ Rab - 1 )))
                if CN_sub > 0.97:
                    CN_SUM = CN_SUM + 1
            else:
                continue
    return CN_SUM

def GCN_multi( atom_id, CN_max = 12 ):   # 计算广义配位数 注意 为了有效处理边界原子情况 计算配位数都用扩胞后的函数计算
    CN_SUM = 0
    atom_id_multi = int(atom_id) + 4 * n
    distances_atom2others = all_distances_multi[atom_id_multi]
    Ra = dict_Pyykko_radii[symbols_multi[atom_id_multi]]
    # if symbols[atom_id] == symbols_multi[atom_id_multi]:
    #     print(symbols[atom_id] + '扩胞正确')
    # else:
    #     print(symbols[atom_id] + '扩胞失败')
    CN_atom_id = []
    CN_atom_id_scaled = []
    for i in range(len(distances_atom2others)):
        if i != int(atom_id_multi):
            Rb = dict_Pyykko_radii[symbols_multi[i]]
            Rab = distances_atom2others[i]
            CN_sub = 1 / (1 + np.exp(-16 * (4/3 * (Ra + Rb)/ Rab - 1 )))
            if CN_sub > 0.97:
                CN_SUM = CN_SUM + 1
                CN_atom_id_scaled.append(i)
                CN_atom_id.append( i - 4 * n )
    # 计算和指定原子的广义配位数
    M_CN_SUM = 0
    GCN = 0
    for j in CN_atom_id_scaled:
        distances_M2others = all_distances_multi[j]
        Ra_M = dict_Pyykko_radii[symbols_multi[j]]
        CN_j = 0
        for k in range(len(distances_M2others)):
            if k != j:
                Rb_M = dict_Pyykko_radii[symbols_multi[k]]
                Rab_M = distances_M2others[k]
                CN_M_sub = 1 / (1 + np.exp(-16 * (4/3 * (Ra_M + Rb_M)/ Rab_M - 1 )))
                if CN_M_sub > 0.97:    # 设置判定成键概率条件
                    CN_j +=1
        GCN += CN_j/CN_max
    print(f"{atom_id} 原子第一配位层配位数为{CN_SUM}，广义配位数为{GCN}")
    return CN_SUM, GCN, CN_atom_id   # 返回指定原子配位数，其配位原子的配位总数以及配位原子列表

# batch_计算相关特征 定义全局变量


global current_path
current_path = os.getcwd()
# features_list = file_name_info_dict[filename]
db_dict = {} # 储存特征的字典
f = open(f'{current_path}/GCN_check_info', 'w')

for foldername, subfolders, filenames in os.walk(current_path):
    with alive_bar( len(subfolders) + 1, title = 'Calculation progress') as bar:
        for i in range(len(subfolders)):
            cwd = os.chdir(subfolders[i])
            # ASE 读取CONTCAR的结构信息以及扩胞操作
            slab = read("CONTCAR")
            symbols = slab.get_chemical_symbols()
            pos = slab.get_positions()
            pos_scaled = slab.get_scaled_positions()
            n = len(slab)
            all_distances = slab.get_all_distances(mic=False)  # 获得所有原子距离向量列表
            an = slab.get_atomic_numbers()  # 获取整数的原子编号数组
            # chemical_formula = slab.get_chemical_formula(mode='hill', empirical=True)  # 基于化学符号将化学配方作为数组
            slab_multi = slab * (3, 3, 1)
            n_multi = len(slab_multi)
            symbols_multi = slab_multi.get_chemical_symbols()
            pos_multi = slab_multi.get_positions()
            all_distances_multi = slab_multi.get_all_distances(mic=False)
            # 计算相关特征
            atom_id = int(sys.argv[1])  # 指定考察广义配位数的原子的id
            CN_max = int(sys.argv[2])  # 指定考察原子的最大配位数
            structure_info = []
            CN_SUM, GCN, CN_atom_id = GCN_multi( atom_id, CN_max = CN_max )

            structure_info.append(CN_SUM)
            structure_info.append(GCN)
            f.write(f'{subfolders[i]} #{str(atom_id)}  原子的配位总数为 {str(CN_SUM)} \n')
            f.write(f'{subfolders[i]} #{str(atom_id)}  原子的配位金属广义配位数为 {str(GCN)} \n')
            f.write(f'{subfolders[i]} #{str(atom_id)}  原子的配位金属列表为 {CN_atom_id} \n')
            bar()
            db_dict.update({subfolders[i]: structure_info})  # update the db_dict
            os.chdir(current_path)
    break

with alive_bar( 2, title='Writing to DATABASE') as bar:
    features_list = ['CN', 'Generalized_CN']
    obj = pd.DataFrame(db_dict, index = features_list ).T  # change the data into DataFrame format
    time.sleep(.001)
    bar()
    obj.to_csv(f'{current_path}/CN_GCN.csv')
    time.sleep(.001)
    bar.text('Task has finished! ')
f.close()

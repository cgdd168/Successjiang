#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description: # usage: python3 batch_CN_DIST.py [file_name] [dopant_id] [O1_id] [02_id] [O3_id]

import os
import copy
from ase import Atoms
from ase.io import read, write
import pandas as pd
import numpy as np
import sys
from alive_progress import alive_bar
import time
from ase.data import covalent_radii
from ase import neighborlist

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

def CN_multi( atom_id ):   # 注意 为了有效处理边界原子情况 计算配位数都用扩胞后的函数计算
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
    # 计算和氧直接相连的金属配位总数
    M_CN_SUM = 0
    for j in CN_atom_id_scaled:
        distances_M2others = all_distances_multi[j]
        M_CN = 0
        Ra_M = dict_Pyykko_radii[symbols_multi[j]]
        for k in range(len(distances_M2others)):
            if k != j:
                Rb_M = dict_Pyykko_radii[symbols_multi[k]]
                Rab_M = distances_M2others[k]
                CN_M_sub = 1 / (1 + np.exp(-16 * (4/3 * (Ra_M + Rb_M)/ Rab_M - 1 )))
                if CN_M_sub > 0.97:
                    M_CN = M_CN + 1
        M_CN_SUM = M_CN_SUM + M_CN
    return CN_SUM, M_CN_SUM, CN_atom_id   # 返回指定原子配位数，其配位原子的配位总数以及配位原子列表

# 计算以掺杂剂为中心的原子距离相邻其他晶胞中心阳离子的距离
def dopant_to_neighbor_cell_center_dist_avg(atom_id): # 输入惨杂元素的初始ASE位置编号
    atom_id_multi = int(atom_id) + 4 * n
    ele = symbols_multi[atom_id_multi]
    distances_atom2others = all_distances_multi[atom_id_multi]
    Ra = dict_Pyykko_radii[ele]
    Rb = dict_Pyykko_radii['V']
    cut_off = 1.5 * (Ra + Rb)     # 设置相邻cell判断截断半径的阈值
    neigbor_len_list = []
    neigbor_id_list = []
    for i in range(len(distances_atom2others)):
        if distances_atom2others[i] < cut_off and symbols_multi[i] == 'V':
            neigbor_len_list.append(distances_atom2others[i])
            neigbor_id_list.append(i - 4 *n)
    center_avg = np.mean(neigbor_len_list)
    return neigbor_id_list, center_avg

# 计算掺杂剂到指定吸附平面的轴向距离，这个特征用来区分表层掺杂或者次表层掺杂
def vertical_dist_dopant_to_ads(atom_id_fixed,atom_id): # atom_id_fixed为定义吸附平面的锚点原子编号；atom_id为输入掺杂原子编号
    """
    描述：point4到point1, point2, point3所在面的距离 point1 为锚点，
    :param point1:数据框的行切片，三维
    :param point4:
    :return:点到面的距离
    """
    atom_id = int(atom_id)
    atom_id_fixed = int(atom_id_fixed)
    pos_fixed = pos[atom_id_fixed]
    pos_fixed_0 = pos_fixed.copy()
    pos_fixed_0[2] +=1.5
    pos_dopant = pos[atom_id]
    n = np.array([0, 0, 1])  # 法向量
    v41 = pos_dopant - pos_fixed_0  # p4到p1的向量
    return abs(v41.dot(n) / n.dot(n))

file_name_info_dict = { 'V2O5_1': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O1', 'SUM_CN_M_O1', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'V2O5_2': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O1', 'SUM_CN_M_O1', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'V2O5_SUB': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O1', 'SUM_CN_M_O1', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'VO2_1': [ 'dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'VO2_2': [ 'dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'VO2_SUB': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O2', 'SUM_CN_M_O2', 'CN_O3', 'SUM_CN_M_O3'],
                        'V2O3_1': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O3', 'SUM_CN_M_O3'],
                        'V2O3_2': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O3', 'SUM_CN_M_O3'],
                        'V2O3_SUB': ['dist_to_cell_center_avg', 'vertical_dist_to_ads_plane', 'CN_O3', 'SUM_CN_M_O3']}
# batch_计算相关特征 定义全局变量
global current_path
current_path = os.getcwd()
filename = sys.argv[1]
features_list = file_name_info_dict[filename]
db_dict = {} # 储存特征的字典
f = open(f'{current_path}/{filename}check_CN_DIST_info', 'w')
count = 0
for foldername, subfolders, filenames in os.walk(current_path):
    for i in range(len(subfolders)):
        count +=1
with alive_bar( 3 * count, title='Calculation progress') as bar:
    for foldername, subfolders, filenames in os.walk(current_path):
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
            chemical_formula = slab.get_chemical_formula(mode='hill')  # 基于化学符号将化学配方作为数组
            slab_multi = slab * (3, 3, 1)
            n_multi = len(slab_multi)
            symbols_multi = slab_multi.get_chemical_symbols()
            pos_multi = slab_multi.get_positions()
            all_distances_multi = slab_multi.get_all_distances(mic=False)
            # 计算相关特征
            dopant_id = int(sys.argv[2])
            atom_id_fixed = int(sys.argv[3])
            O_id_list = [int(i) for i in sys.argv[3:]]
            structure_info = []
            neigbor_id_list, center_avg = dopant_to_neighbor_cell_center_dist_avg(dopant_id)  # 计算以掺杂剂为中心的原子距离相邻其他晶胞中心阳离子的距离
            structure_info.append(center_avg)
            f.write(f'{subfolders[i]} #{str(dopant_id)}  原子的邻近晶胞列表为 {str(neigbor_id_list)} \n')
            f.write(f'{subfolders[i]} #{str(dopant_id)}  原子的邻近晶胞平均距离为 {str(center_avg)} \n')
            bar()
            vertical_dist_to_ads = vertical_dist_dopant_to_ads(atom_id_fixed, dopant_id) # 计算掺杂剂到指定吸附平面的轴向距离，这个特征用来区分表层掺杂或者次表层掺杂
            structure_info.append(vertical_dist_to_ads)
            f.write(f'{subfolders[i]} #{str(dopant_id)}  掺杂剂到指定吸附平面的轴向距离 {str(vertical_dist_to_ads)} \n')
            bar()
            for O_id in O_id_list:
                O_CN_SUM, O_M_CN_SUM, O_CN_atom_id_list = CN_multi(O_id)
                structure_info.append(O_CN_SUM)
                structure_info.append(O_M_CN_SUM)
                f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位总数为 {str(O_CN_SUM)} \n')
                f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子配位总数为 {str(O_M_CN_SUM)} \n')
                f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子列表为 {O_CN_atom_id_list} \n')
            bar()
            db_dict.update({subfolders[i]: structure_info})  # update the db_dict
            os.chdir(current_path)
with alive_bar( 2, title='Writing to DATABASE') as bar:
    obj = pd.DataFrame(db_dict, index = features_list ).T  # change the data into DataFrame format
    time.sleep(.001)
    bar()
    obj.to_csv(f'{current_path}/{filename}_CN_DIST.csv')
    time.sleep(.001)
    bar(text='Task has finished! ')
f.close()

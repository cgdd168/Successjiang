#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description:
import numpy as np
from sklearn import preprocessing
data = np.loadtxt('TRAIN1.txt', delimiter='\t')
z_scaler= preprocessing.StandardScaler() # 建立 StandardScaler 对象
z_data = z_scaler.fit_transform(data) # 用 StandardScaler 对象对数据进行标准化处理
np.savetxt('data_z_score.txt',z_data)
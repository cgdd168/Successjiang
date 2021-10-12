#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description:
import pandas as pd
import numpy as np
import os
import seaborn as sns
from alive_progress import alive_bar
import matplotlib.pyplot as plt

global current_path
current_path = os.getcwd()
db_dict = {} # 储存结果的字典
data = {}
files = os.listdir(current_path) # 遍历filepath下所有文件，包括子目录
i=0
with alive_bar(28, title='Total progress of the task') as bar:
    for rung in files:
        fi_d = os.path.join(current_path, rung)
        if os.path.isdir(fi_d):  # isdir和isfile参数必须跟绝对路径
            os.chdir(fi_d)
            Dimensions = os.listdir(fi_d)
            for D in Dimensions:
                D_path = os.path.join(fi_d, D)
                if os.path.isdir(D_path):
                    os.chdir(D_path)
                    print(f"当前目录是{D_path}")
                    os.system('grep "Total RMSE,MaxAE" SISSO.out >result')
                    with open("result",'r') as output:
                        results=output.readlines()
                        RMSE = float(results[-1].strip().split(":")[1].strip().split()[0])
                        bar()
                        db_dict.update({D: RMSE})
                        Ser = pd.Series(db_dict)
                        Ser.index.name = 'Dimension'
                        Ser.name = f'{rung}'
            data.update({rung: Ser})

frame = pd.DataFrame(data)
os.chdir(current_path)
frame.to_csv(f'{current_path}/Total_Training_Results.csv')

# 做出统计分布图
sns.set_style('ticks')
sns.set_palette('Set1')
fig = frame.plot(figsize=(6,4)) # figsize：创建图表窗口，设置窗口大小
plt.title('Total_train_Set_Results')  # 图名
plt.xlabel('Descriptor Dimensionn')  # x轴标签
plt.xlim(0,11)
plt.ylabel('RMSE(eV)') # y轴标签
plt.savefig('Training_Total_Results.svg', bbox_inches='tight')








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
import sys
import matplotlib.pyplot as plt
from alive_progress import alive_bar


def RMSE_MaxAE():
    files_total = os.listdir(current_path)
    if 'CV_RMSE' in files_total:
        with open('CV_RMSE', 'r') as CV_RMSE:
            n = CV_RMSE.read()
            if n == 0:
                os.system('SISSO_predict ')
                os.system('grep RMSE predict_Y.out> CV_RMSE')
                with open('CV_RMSE', 'r') as CV:
                    results = CV.readlines()
                    output=[]
                    RMSE = float(results[-1].strip().split(":")[1].strip().split()[0])
                    MaxAE = float(results[-1].strip().split(":")[1].strip().split()[1])
                    output.append(RMSE)
                    output.append(MaxAE)
                    return output
    else:
        os.system('SISSO_predict')
        os.system('grep RMSE predict_Y.out> CV_RMSE')
        with open('CV_RMSE', 'r') as CV:
            results = CV.readlines()
            output = []
            RMSE = float(results[-1].strip().split(":")[1].strip().split()[0])
            MaxAE = float(results[-1].strip().split(":")[1].strip().split()[1])
            output.append(RMSE)
            output.append(MaxAE)
    return output

global current_path
current_path = os.getcwd()
db_dict = {} # 储存结果的字典
files = os.listdir(current_path) # 遍历filepath下所有文件，包括子目录
i=0
with alive_bar(10, title='Total progress of the task') as bar:
    for fi in files:
        fi_d = os.path.join(current_path, fi)
        if os.path.isdir(fi_d):  # isdir和isfile参数必须跟绝对路径
            i+=1
            os.chdir(fi_d)
            print(f"正在打开{fi_d}文件夹")
            results = RMSE_MaxAE()
            bar()
            db_dict.update({i: results})  # update the dict
            os.chdir(current_path)

# 绘制结果分析图
# font_options = {'family' : 'monospace',
#                 'weight' : 'bold',
#                 'size': 'small'}
# plt.rc('font', **font_options)

sns.set_style('ticks')
sns.set_palette('Set1')
obj = pd.DataFrame(db_dict, index = ['RMSE', 'MaxAE'], columns=np.arange(1,11,1)).T
obj=obj.sort_values('RMSE')
obj.index=range(1,11)

RMSE_mean = round(obj.mean()[0], 2)
MaxAE_mean = round(obj.mean()[1], 2)
with open('RMSE_MaxAE_avg','w') as output:
    output.write(f'RMSE_mean is {RMSE_mean}\n')
    output.write(f'MaxAE_mean is {MaxAE_mean}\n')
fig = obj.plot(figsize=(6,4)) # figsize：创建图表窗口，设置窗口大小
plt.title('CV_10_Results')  # 图名
plt.xlabel('Iteration')  # x轴标签
plt.xlim(0,11)
plt.xticks(np.linspace(1,10,10, endpoint= True))
plt.ylabel('RMSE/MaxAE(eV)') # y轴标签
plt.axhline(y=RMSE_mean,color="red",ls="--") # 在图上添加RMSE均值线
plt.text(0,RMSE_mean+0.05, f'RMSE-avg={RMSE_mean}')
plt.axhline(y=MaxAE_mean,color="blue",ls="--") # 在图上添加MaxAE均值线
plt.text(0,MaxAE_mean+0.05, f'MaxAE-avg={MaxAE_mean}')
plt.savefig('CV_10_Results.svg', bbox_inches='tight')
obj.to_csv(f'{current_path}/CV_10_Results.csv')


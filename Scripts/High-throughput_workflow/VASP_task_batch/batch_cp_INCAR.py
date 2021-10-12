#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description:
import os
import shutil
count = 0
global current_path
current_path = os.getcwd()
for foldername, subfolders, filenames in os.walk(current_path):
    for i in range(len(subfolders)):
        os.chdir(subfolders[i])
        os.system('cp ../INCAR .')
        os.system('cp ../KPOINTS .')
        os.system('cp ../script .')
        os.system('cif2pos.py *.cif')
        os.system('mkpotcar')
        print('continued files required of ' + subfolders[i] + ' has been prepared!')
        os.chdir(current_path)
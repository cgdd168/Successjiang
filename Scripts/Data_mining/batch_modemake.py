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
from alive_progress import alive_bar
from ase.io import read, write

global current_path
# global Filename
current_path = os.getcwd()

db_dict = {}
for foldername, subfolders, filenames in os.walk(current_path):
    with alive_bar(len(subfolders), title='Total progress of the task') as bar:
        for i in range(len(subfolders)):
            cwd = os.chdir(subfolders[i])

            slab = read("CONTCAR")
            n = len(slab)
            os.system(f'grep "33 f/i" OUTCAR -A {n + 1} >freq33')
            os.system('modemake.py freq33 0.1')
            os.system('mkdir dim')
            os.system('cp ~/INCAR_template/dim/INCAR ./dim')
            os.system('cp ./KPOINTS ./dim')
            os.system('cp ./POTCAR ./dim')
            os.system('cp ./POSCAR ./dim')
            os.system('cp ~/INCAR_template/dim/script ./dim')
            os.system('cp ./MODECAR ./dim')
            bar()
            os.chdir(current_path)
        break

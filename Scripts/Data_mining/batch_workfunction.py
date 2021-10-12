#! python3
# batch_DOS.py - copy INCAR, KPOINTS, script in the current working folder into the destination folder.
#                copy CONTCAR to POSCAR of every destination folder, respectively.
# usage:  prepare INCAR, KPOINTS, script files for DOS and all the folders with structures optimized (i.e. CONTCAR) in the current working folder.
#         run python3 batch_DOS.py
# author: Jiangcg

import os
import shutil

count = 0
global current_path
current_path = os.getcwd()
for foldername, subfolders, filenames in os.walk(current_path):
    for i in range(len(subfolders)):
        os.chdir(subfolders[i])
        sub_file = os.getcwd()
        shutil.copy(current_path + '/INCAR', sub_file)
        shutil.copy(current_path + '/KPOINTS', sub_file)
        os.system(f'cp {sub_file}/CONTCAR {sub_file}/POSCAR ; rm REPORT CHG* DOSCAR EIGENVAL IBZKPT PCDAT PROCAR WAVECAR XDATCAR vasprun.xml FORCECAR NODEIB script')
        print('work-function files required of ' + subfolders[i] + ' has been created!')
        count += 1
        os.chdir(current_path)

print('ALL ' + str(count) + ' work function FILES HAVE DONE! GOOD LUCK!')

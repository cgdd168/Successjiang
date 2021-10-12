element = ['Ag','Au','Cd','Ce','Co','Cr','Cu','Fe','Hf','Hg','Ir','Mn','Mo','Nb','Ni','Os','Pd','Pt','Re','Rh','Ru','Sc','Ta','Tc','Ti','W','Y','Zn','Zr']

import os
import shutil

current_path = os.getcwd()
for i, n in enumerate(element) :

    ele_file = current_path + '/' + str(n)
    shutil.copy(current_path + '/INCAR', ele_file)
    shutil.copy(current_path + '/KPOINTS', ele_file)
    shutil.copy(current_path + '/script', ele_file)
    os.system(f'cp {ele_file}/CONTCAR {ele_file}/POSCAR')
    os.system(f'rm {ele_file}/OUTCAR')

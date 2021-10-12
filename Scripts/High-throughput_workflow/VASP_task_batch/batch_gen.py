element = ['Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
'Ce', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']

import os
import shutil

current_path = os.getcwd()

for ele in element:
    old_str = 'Lu'
    new_str = ele

    ele_file = current_path + '/' + str(new_str)
    term = ele_file + '/POSCAR'
    
    if not os.path.exists(ele_file):
        os.mkdir(ele_file)

    file_data = ""
    with open(current_path + "/POSCAR", "r") as f:
        for line in f:
            line = line.replace(old_str,new_str)
            file_data += line
    
    with open(term, "w") as f:
        f.write(file_data)
    
    shutil.copy(current_path + '/INCAR', ele_file)
    shutil.copy(current_path + '/KPOINTS', ele_file)
    shutil.copy(current_path + '/script', ele_file)
    os.chdir(ele_file)
    os.system('mkpotcar')

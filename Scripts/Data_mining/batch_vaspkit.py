#! /usr/bin/python3
# batch_vaspkit.py - do vaspkit analysis from all subfolders in the current working folder.
# written by - Jiangcg & Wusc
# usage: python3 batch_vaspkit.py i
#        what's the i meaning?
#        1 - D band center calculation for all atoms
#        the d-band center for every atoms and the total d-band center is writen in D_BAND_CENTER file
#
#        2 - Get the total DOS Input:
#        Input: None.
#        Output:
#        TDOS.dat contains the total DOS;
#        ITDOS.dat contains the integral total DOS.
#        Spin up and down are written in one file.
#        Example: Total dos of Θ-Al2O3 unit cell
#
#        3 - Output projected DOS for every elements to separate files.
#        Input: None.
#        Output: PDOS_Elements_UP (_DW).dat, pdos file for each elements (sum of all atoms of each element.)
#        IPDOS_Elements_UP (_DW).datSpin up and down are writen in different files.
#

import os, sys

vaspkit_function = {'1': ['d-band center of', '503\n', 'n\n', 'D_BAND_CENTER'], '2': ['Total DOS', '111\n', '0\n', 'TDOS.dat'],
                    '3': ['projected DOS for every elements','113\n','0\n', 'PDOS_Element_UP (_DW).dat']}
count = 0
global current_path
current_path = os.getcwd()
for foldername, subfolders, filenames in os.walk(current_path):
    for i in range(len(subfolders)):
        cwd = os.chdir(subfolders[i])
        if len(sys.argv) == 2 and int(sys.argv[1]) >=1 and  int(sys.argv[1]) <= 3:
            a = sys.argv[1]
            os.system("echo -e '%s%s'|vaspkit > /dev/null&" % (vaspkit_function[a][1], vaspkit_function[a][2]))
            if sys.argv[1] != '3' and vaspkit_function[a][3] in os.listdir(cwd):
                print( vaspkit_function[a][0] + ' ' + subfolders[i] + ' has been calculated! ')
                count += 1
            elif sys.argv[1] == '3':
                print(vaspkit_function[a][0] + ' ' + subfolders[i] + ' has been calculated! ')
                count += 1
            else:
                print( vaspkit_function[a][0] + ' ' + subfolders[i] + ' has NOT done !!!! ATTENTION!')
        os.chdir(current_path)

print( str(count) + ' folders ' + vaspkit_function[sys.argv[1]][0] + ' have been done! Good luck!' )
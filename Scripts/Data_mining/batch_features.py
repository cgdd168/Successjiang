#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File Name: 
# Version: V1.0
# @Author: Jiangcg
# Email: successjiang@tju.edu.cn
# Organization: Tianjin University
# @Date:
# Description:
'''
This script has automatic batch analysis functions from output files after VASP DOS calculations:
Electronic structure information:
 1) Workfunction calculation
 2) Bader charge analysis and calculate the net charge for specific atom (ASE_id)
 3) Band-center (s,p,d,f band center for specific atom)
 4) Band-width (s,p,d,f band center for specific atom)
 5) Unoccupied band center (p,d,f band center above the fermi level for specific atom)
 6) Band-filling

Geometric structure information:
 1) Coordination number calculation
 2) GCN (not including in this script, but can be obtained in another script in the same folder)

Please Contact successjiang@tju.edu.cn if you have questions.

Note1: Please make sure that you have installed VASPKIT software program correctly, because this is used to obtained
vacuum level when calculating Work Function.

Note2: Coordination number calculation bases on the experimental calculation equation from original paper about DFT-D3.
  This means
  1) you have to test the P_band number (0~1) which specifies the covalent bonding probability between two atoms.
  2) Recommended P_band is: 0.97 for transition metal oxide system; 0.80 for transition metal system.
  3) All atom geometric information read from CONTCAR using ASE Module.
  4) In order to deal with the periodic boundary atoms accurately, a 3*3 extending supercell is used atomatically.

Note3: All analysis outputs are written in the main folder in the form of CSV using Pandas Module.

# Usage: python3 batch_features.py [O_ASE_id] [M1_ASE_id] [M2_ASE_id] ...

'''
import pandas as pd
import numpy as np
import os
import sys
import bisect
from ase.io import read, write
from scipy.integrate import simps
from alive_progress import alive_bar


dict_Pyykko_radii = {'O': 0.63, 'Sc': 1.49, 'Ti': 1.36, 'V': 1.34, 'Cr': 1.22, 'Mn': 1.19,
                     'Fe': 1.16, 'Co': 1.11, 'Ni': 1.1, 'Cu': 1.12, 'Zn': 1.18, 'Y': 1.63,
                     'Zr': 1.54, 'Nb': 1.47, 'Mo': 1.38, 'Tc': 1.28, 'Ru': 1.25, 'Rh': 1.25,
                     'Pd': 1.2, 'Ag': 1.28, 'Cd': 1.36, 'Ce': 1.63, 'Hf': 1.52, 'Ta': 1.46,
                     'W': 1.37, 'Re': 1.31, 'Os': 1.29, 'Ir': 1.22, 'Pt': 1.23, 'Au': 1.24,
                     'Hg': 1.33, 'La': 1.80, 'H': 0.32, 'Li': 1.33, 'Be': 1.02, 'B': 0.85,
                     'C': 0.75, 'N': 0.71, 'k': 1.96, 'F': 0.64, 'Ne': 0.67, 'Na': 1.55,
                     'Mg': 1.39, 'Al': 1.26, 'Si': 116, 'P': 1.11, 'S':1.03, 'Cl': 0.99,
                     'Ca': 1.71, 'Ga': 1.24, 'Ge': 1.21, 'As': 1.21, 'Se': 1.16, 'Br': 1.14,
                     'Kr': 1.17, 'Rb': 2.10, 'Sr': 1.85, 'In': 1.42, 'Sn': 1.40, 'Te': 1.36,
                     'I': 1.33, 'Xe': 1.31, 'Cs': 2.32, 'Ba': 1.96, 'Tl': 1.44, 'Pb': 1.44,
                     'Bi': 1.51, 'Po': 1.45, 'At': 1.47, 'Rn': 1.42, 'Fr': 2.23, 'Ra': 2.01,
                     'Rf': 1.57, 'Db': 1.49, 'Sg': 1.43, 'Bh': 1.41, 'Hs': 1.34, 'Mt': 1.29,
                     'Ds': 1.28, 'Rg': 1.21, 'Pr': 1.76, 'Nd': 1.74, 'Pm': 1.73, 'Sm': 1.72,
                     'Eu': 1.68, 'Gd': 1.69, 'Tb': 1.68, 'Dy': 1.67, 'Ho': 1.66, 'Er': 1.65,
                     'Tm': 1.64, 'Yb': 1.70, 'Lu': 1.62, 'Ac': 1.86, 'Th': 1.75, 'Pa': 1.69,
                     'U': 1.70, 'Np': 1.71, 'Pu': 1.72, 'Am': 1.66, 'Cm': 1.66, 'Bk':1.68,
                     'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67, 'Md': 1.73,'No': 1.76, 'Lr': 1.61 }

P_bond = 0.97 # 成键概率判定，需要根据具体体系测试，金属氧化物可以按0.97开始；金属按0.80开始

# 定义功函数计算函数
def work_function():
    Abs_path = os.path.abspath('vacumm_fermi_level')
    with open(Abs_path,'r') as inputFiles:
        InputFiles = inputFiles.readlines()
        vacumm_level = float(InputFiles[1].split(':')[1].split('\n')[0])
        fermi_level = float(InputFiles[2].split(':')[1].split('X')[0])
        WF = vacumm_level - fermi_level
        print(Abs_path, ' work function', ' is', str(WF))
    return WF

# 定义bader电荷-净电荷计算函数
def bader_charge(atom_id):   # atom_id = atom_id_in ase + 1
#            try:
#               if 'ACF.dat' in os.listdir('./'):

    os.system("sed -n '%s,%sp' ACF.dat > Charge%s " % (atom_id+2, atom_id+2, atom_id) )
    os.system("sed -n '6,7p' POSCAR >> Charge%s " % atom_id )
    os.system("grep -w ZVAL POTCAR >> Charge%s " % atom_id)
    Abs_path = os.path.abspath(f'Charge{atom_id}')
    with open(Abs_path, 'r') as inputFiles: # InputFiles = open('ACF.dat', 'r').readlines()
        InputFiles = inputFiles.readlines()
        ZVAL_keys = []
        ZVAL_values = []
        POSCAR_atom_ID = InputFiles[2].strip().split()
        n = 1
        for i in range(len(POSCAR_atom_ID)):   # 构建一个包含价电子的字典，以POSCAR中元素编号为键
            n = n + int(POSCAR_atom_ID[i])
            ZVAL_keys.append(n)
            ZVAL_values.append(float(InputFiles[i+3].split('L   =')[1].split('mass')[0]))
        ZVAL_dict = {key: value for key, value in zip(ZVAL_keys,ZVAL_values )}
        Charge = float(InputFiles[0].strip().split()[4]) # 提取ACF.dat中特定原子的bader电荷
        pos_0 = bisect.bisect(list(ZVAL_dict.keys()), atom_id)    # 通过二分法找到特定原子在POSCAR中的对应元素位置
        pos_1 = list(ZVAL_dict.keys())[pos_0]                     # 通过二分法找到特定原子在POSCAR中的对应元素位置
        atom_zval = ZVAL_dict[pos_1]
        Bader_Charge = atom_zval - Charge
    return Bader_Charge

# @jit(nopython=True)
def pdos_column_names(lmax, ispin):  # 首次调用时，函数被编译为机器代码
    if lmax == 2:
        names = ['s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2']
        # names = [ 's', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2' ]
    elif lmax == 3:
        names = ['s', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2',
                 'f_y(3x2-y2)', 'f_xyz', 'f_yz2', 'f_z3', 'f_xz2', 'f_z(x2-y2)', 'f_x(x2-3y2)']
        # names = [ 's', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2','f-3',
        #         'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']
    else:
        raise ValueError('lmax value not supported')

    if ispin == 2:
        all_names = []
        for n in names:
            all_names.extend(['{}_up'.format(n), '{}_down'.format(n)])
    else:
        all_names = names
    all_names.insert(0, 'energy')
    return all_names


class Doscar:
    '''
    Contains all the data in a VASP DOSCAR file, and methods for manipulating this.
    '''

    number_of_header_lines = 6

    def __init__(self, filename, ispin=2, lmax=2, lorbit=11, spin_orbit_coupling=False, read_pdos=True, species=None):
        '''
        Create a Doscar object from a VASP DOSCAR file.
        Args:
            filename (str): Filename of the VASP DOSCAR file to read.
            ispin (optional:int): ISPIN flag.
                Set to 1 for non-spin-polarised or 2 for spin-polarised calculations.
                Default = 2.
            lmax (optional:int): Maximum l angular momentum. (d=2, f=3). Default = 2.
            lorbit (optional:int): The VASP LORBIT flag. (Default=11).
            spin_orbit_coupling (optional:bool): Spin-orbit coupling (Default=False).
            read_pdos (optional:bool): Set to True to read the atom-projected density of states (Default=True).
            species (optional:list(str)): List of atomic species strings, e.g. [ 'Fe', 'Fe', 'O', 'O', 'O' ].
                Default=None.
        '''
        self.filename = filename
        self.ispin = ispin
        self.lmax = lmax
        self.spin_orbit_coupling = spin_orbit_coupling
        if self.spin_orbit_coupling:
            raise NotImplementedError('Spin-orbit coupling is not yet implemented')
        self.lorbit = lorbit
        self.pdos = None
        self.species = species
        self.read_header()
        self.read_total_dos()
        if read_pdos:
            try:
                self.read_projected_dos()
            except:
                raise
        # if species is set, should check that this is consistent with the number of entries in the
        # projected_dos dataset

    # @jit(nopython=True)
    @property
    def number_of_channels(self):
        if self.lorbit == 11:
            return {2: 9, 3: 16}[self.lmax]
        raise NotImplementedError

    # @jit(nopython=True)
    def read_header(self):
        self.header = []
        with open(self.filename, 'r') as file_in:
            for i in range(Doscar.number_of_header_lines):
                self.header.append(file_in.readline())
        self.process_header()

    # @jit(nopython=True)
    def process_header(self):
        self.number_of_atoms = int(self.header[0].split()[0])
        self.number_of_data_points = int(self.header[5].split()[2])
        self.efermi = float(self.header[5].split()[3])

    # @jit(nopython=True)
    def read_total_dos(self):  # assumes spin_polarised
        start_to_read = Doscar.number_of_header_lines
        df = pd.read_csv(self.filename,
                         skiprows=start_to_read,
                         nrows=self.number_of_data_points,
                         delim_whitespace=True,
                         names=['energy', 'up', 'down', 'int_up', 'int_down'],
                         index_col=False)
        self.energy = df.energy.values
        df.drop('energy', axis=1)
        self.tdos = df

    # @jit(nopython=True)
    def read_atomic_dos_as_df(self, atom_number):  # currently assume spin-polarised, no-SO-coupling, no f-states
        assert atom_number > 0 & atom_number <= self.number_of_atoms
        start_to_read = Doscar.number_of_header_lines + atom_number * (self.number_of_data_points + 1)
        df = pd.read_csv(self.filename,
                         skiprows=start_to_read,
                         nrows=self.number_of_data_points,
                         delim_whitespace=True,
                         names=pdos_column_names(lmax=self.lmax, ispin=self.ispin),
                         index_col=False)
        return df.drop('energy', axis=1)

    # @jit(nopython=True)
    def read_projected_dos(self):
        """
        Read the projected density of states data into """
        pdos_list = []
        for i in range(self.number_of_atoms):
            df = self.read_atomic_dos_as_df(i + 1)
            pdos_list.append(df)
        # self.pdos  =   pdos_list
        self.pdos = np.vstack([np.array(df) for df in pdos_list]).reshape(
            self.number_of_atoms, self.number_of_data_points, self.number_of_channels, self.ispin)

    # @jit(nopython=True)
    def pdos_select(self, atoms=None, spin=None, l=None, m=None):
        """
        Returns a subset of the projected density of states array.
        """
        valid_m_values = {'s': [],
                          'p': ['x', 'y', 'z'],
                          'd': ['xy', 'yz', 'z2-r2', 'xz', 'x2-y2'],
                          'f': ['y(3x2-y2)', 'xyz', 'yz2', 'z3', 'xz2', 'z(x2-y2)', 'x(x2-3y2)']}
        if not atoms:
            atom_idx = list(range(self.number_of_atoms))
        else:
            atom_idx = atoms
        to_return = self.pdos[atom_idx, :, :, :]
        if not spin:
            spin_idx = list(range(self.ispin))
        elif spin is 'up':
            spin_idx = [0]
        elif spin is 'down':
            spin_idx = [1]
        elif spin is 'both':
            spin_idx = [0, 1]
        else:
            raise ValueError
        to_return = to_return[:, :, :, spin_idx]

        if not l:
            channel_idx = list(range(self.number_of_channels))
        elif l == 's':
            channel_idx = [0]
        elif l == 'p':
            if not m:
                channel_idx = [1, 2, 3]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['p']) if v in m]
        elif l == 'd':
            if not m:
                channel_idx = [4, 5, 6, 7, 8]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['d']) if v in m]
        elif l == 'f':
            if not m:
                channel_idx = [9, 10, 11, 12, 13, 14, 15]
            else:
                channel_idx = [i for i, v in enumerate(valid_m_values['f']) if v in m]
        else:
            raise ValueError

        return to_return[:, :, channel_idx, :]

    def pdos_sum(self, atoms=None, spin=None, l=None, m=None):
        return np.sum(self.pdos_select(atoms=atoms, spin=spin, l=l, m=m), axis=(0, 2, 3))

# 计算指定原子的带中心
def band_center_width(atom, orb, subfolder):
    # calculation of band center
    # Open doscar
    dosfile = 'DOSCAR'
    doscar = Doscar(dosfile, ispin=2, lmax=3, lorbit=11)  # calculation setting

    # erange_down_limit = float(doscar.header[5].split()[1])
    # erange_up_limit = float(doscar.header[5].split()[0])
    # Set atoms for integration
    a = int(atom)
    atom = [a]  # calculated atom ordinal
    orb = orb
    up = doscar.pdos_sum(atom, spin='up', l=orb)
    down = doscar.pdos_sum(atom, spin='down', l=orb)
    both = doscar.pdos_sum(atom, spin='both', l=orb)

    # Set intergrating range
    efermi = doscar.efermi - doscar.efermi
    energies = doscar.energy - doscar.efermi
    emin, emax = energies[0], energies[-1]
    erange = (emin, emax)  # integral energy range
    # erange = (efermi-21, efermi+6)      # integral energy range
    emask = (energies <= erange[-1])

    # Calculating center of the orbital specified above in line 184
    x = energies[emask]
    y1 = up[emask]
    y2 = down[emask]
    y3 = both[emask]

    dbc_up = simps(y1 * x, x) / simps(y1, x)
    dbc_down = simps(y2 * x, x) / simps(y2, x)
    dbc_avg = (dbc_up + dbc_down)/2
    # dbc_avg = simps(y3 * x, x) / simps(y3, x)
    band_width = np.sqrt( simps(y3 * x * x, x) / simps(y3, x) )
    print(subfolder, ' ', atom, f' atom {orb} band center[up, down, 0.5*(up_down)] and band_width] is ', dbc_up, dbc_down, dbc_avg, band_width)
    return dbc_avg, band_width

# 计算指定原子的空带中心
def unoccupied_band_center(atom, orb, subfolder):
    # calculation of band center
    # Open doscar
    dosfile = 'DOSCAR'
    doscar = Doscar(dosfile, ispin=2, lmax=3, lorbit=11)  # calculation setting

    # erange_down_limit = float(doscar.header[5].split()[1])
    # erange_up_limit = float(doscar.header[5].split()[0])
    # Set atoms for integration
    a = int(atom)
    atom = [a]  # calculated atom ordinal
    orb = orb
    up = doscar.pdos_sum(atom, spin='up', l=orb)
    down = doscar.pdos_sum(atom, spin='down', l=orb)
    both = doscar.pdos_sum(atom, spin='both', l=orb)

    # Set intergrating range
    efermi = doscar.efermi - doscar.efermi
    energies = doscar.energy - doscar.efermi
    emin, emax = energies[0], energies[-1]
    erange = (efermi, emax)  # integral energy range
    # erange = (efermi-21, efermi+6)      # integral energy range
    emask = (energies >= erange[0])

    # Calculating center of the orbital specified above in line 184
    x = energies[emask]
    y1 = up[emask]
    y2 = down[emask]
    y3 = both[emask]
    dbc_up = simps(y1 * x, x) / simps(y1, x)
    dbc_down = simps(y2 * x, x) / simps(y2, x)
    d_un_center = (dbc_up + dbc_down)/2
    # d_un_center = simps(y3 * x, x) / simps(y3, x)
    # band_width = np.sqrt( simps(y3 * x * x, x) / simps(y3, x) )
    # dbc = []
    # dbc.append(dbc_avg)
    print(subfolder, ' ', atom, f' atom fermi_above d-band center[avg] is ', d_un_center)
    return d_un_center

# 计算指定原子的带填充程度(%)-电子能量为参考
def band_filling(atom, orb, subfolder):
    # calculation of band center
    # Open doscar
    dosfile = 'DOSCAR'
    doscar = Doscar(dosfile, ispin=2, lmax=3, lorbit=11)  # calculation setting

    # erange_down_limit = float(doscar.header[5].split()[1])
    # erange_up_limit = float(doscar.header[5].split()[0])
    # Set atoms for integration
    a = int(atom)
    atom = [a]  # calculated atom ordinal
    orb = orb
    up = doscar.pdos_sum(atom, spin='up', l=orb)
    down = doscar.pdos_sum(atom, spin='down', l=orb)
    both = doscar.pdos_sum(atom, spin='both', l=orb)

    # Set intergrating range
    efermi = doscar.efermi - doscar.efermi
    energies = doscar.energy - doscar.efermi
    emin, emax = energies[0], energies[-1]
    erange1 = (emin, emax)  # integral energy range
    erange2 = (emin, efermi)  # integral energy range
    # erange = (efermi-21, efermi+6)      # integral energy range
    # emask1 = (energies >= erange1[0])
    emask1 = (energies <= erange1[-1])
    # emask2 = (energies >= erange2[0])
    emask2 = (energies <= erange2[-1])

    # Calculating center of the orbital specified above in line 184
    x1 = energies[emask1]
    x2 = energies[emask2]

    y3_1 = both[emask1]
    y3_2 = both[emask2]

    Band_filling = simps(y3_2 , x2) / simps(y3_1 , x1) # 费米能级以下电子态的积分占全部的百分比
    print(subfolder, ' ', atom, f' atom {orb} band filling is ', Band_filling)
    # band_center_log.write(str(subfolder + ' ' + atom + ' atom fermi_above d-band center[avg] is ' + dbc_avg))
    return Band_filling
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
                if CN_sub > P_bond:
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
            if CN_sub > P_bond:
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
                if CN_M_sub > P_bond:
                    M_CN = M_CN + 1
        M_CN_SUM = M_CN_SUM + M_CN
    return CN_SUM, M_CN_SUM, CN_atom_id   # 返回指定原子配位数，其配位原子的配位总数以及配位原子列表


# 执行计算任务的代码
global current_path
current_path = os.getcwd()
db_dict = {}
filename = sys.argv[1]
f = open(f'{current_path}/NEW_BAND_CENTER_{filename}_features_check_info', 'w')
for foldername, subfolders, filenames in os.walk(current_path):
    with alive_bar(13 * len(subfolders), title='Total progress of the task') as bar:
        for i in range(len(subfolders)):
            cwd = os.chdir(subfolders[i])
            OUTPUT=[]                                     # 输出列表
            os.system('dobader')                          # 执行dobader脚本
            bar()

            os.system(  "echo -e '426\n3\n'|vaspkit | grep Vacuum-Level > vacumm_fermi_level ; grep E-fermi OUTCAR >> vacumm_fermi_level" )
            bar()
            # 计算吸附位点氧的bader电荷
            O_CN_number = len(sys.argv)-2   # 定义氧离子的配位数
            BADER_CHARGE_O = bader_charge(atom_id= int(sys.argv[1]))  # 计算氧原子的bader电荷
            OUTPUT.append(BADER_CHARGE_O)
            bar()
            # 计算吸附位点氧的直接配位的金属阳离子的平均bader电荷
            BADER_CHARGE_M = 0
            for j in range(O_CN_number):
                BADER_CHARGE_M += bader_charge(atom_id = int(sys.argv[j+2]))  # 计算周围配位金属原子的bader电荷
            BADER_CHARGE_M_avg = BADER_CHARGE_M/O_CN_number
            OUTPUT.append(BADER_CHARGE_M_avg)
            bar()
            # 计算吸附位点O的p带中心以及p带填充
            p_band_center, p_band_width = band_center_width(atom=int(sys.argv[1]), orb='p', subfolder=subfolders[i])
            p_band_filling = band_filling(atom=int(sys.argv[1]), orb='p', subfolder=subfolders[i])
            OUTPUT.append(p_band_center)
            OUTPUT.append(p_band_filling)
            bar()

            # 计算吸附位点O周围金属阳离子的平均d带中心以及平均d带宽度
            D_band_center = []
            D_band_width = []
            D_un_center = []
            D_filling = []
            for k in range(O_CN_number):
                d_band_center, d_band_width = band_center_width(atom=int(sys.argv[k+2]), orb='d', subfolder=subfolders[i])
                d_un_center = unoccupied_band_center(atom=int(sys.argv[k+2]), orb='d', subfolder=subfolders[i])
                D_band_filling = band_filling(atom=int(sys.argv[k+2]), orb='d', subfolder=subfolders[i])

                D_un_center.append(d_un_center)
                D_band_center.append(d_band_center)
                D_band_width.append(d_band_width)
                D_filling.append(D_band_filling)

            D_band_center_avg = np.mean(D_band_center)
            D_band_width_avg = np.mean(D_band_width)
            D_un_center_avg = np.mean(D_un_center)
            D_un_center_min = np.min(D_un_center)
            D_band_filling_avg = np.mean(D_filling)

            OUTPUT.append(D_band_center_avg)
            bar()

            OUTPUT.append(D_band_width_avg)
            bar()

            OUTPUT.append(D_un_center_avg)
            bar()

            OUTPUT.append(D_un_center_min)
            bar()

            OUTPUT.append(D_band_filling_avg)
            bar()

            # 计算当前体系表面功函数
            WF = work_function()               # 计算当前体系表面功函数
            OUTPUT.append(WF)
            bar()

            # 计算金属阳离子配位总数以及氧的配位总数

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

            # 计算氧配位总数以及金属原子配位总数
            O_CN_SUM, O_M_CN_SUM, O_CN_atom_id_list = CN_multi(int(sys.argv[1]))

            OUTPUT.append(O_M_CN_SUM)
            bar()

            OUTPUT.append(O_CN_SUM)
            bar()
            O_id = sys.argv[1]
            # 写计算特征细节信息以供校验计算准确性
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位总数为 {str(O_CN_SUM)} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子配位总数为 {str(O_M_CN_SUM)} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子计算列表为 {O_CN_atom_id_list} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子指定列表为 {sys.argv[2:]} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子d带中心分别为为 {D_band_center} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子d带宽度分别为为 {D_band_width} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子d空带中心分别为为 {D_un_center} \n')
            f.write(f'{subfolders[i]} #{str(O_id)}  氧的配位金属阳离子d带填充度别为为 {D_filling} \n')
            db_dict.update({subfolders[i]: OUTPUT})  # update the dict
            os.chdir(current_path)
    break
f.close()

with alive_bar( 2, title='Writing to DATABASE') as bar:

    features_list = ['Q_O', 'Q_M', 'p-center', 'p-filling', 'd-center', 'd-width', 'd-un-center', 'd-un-center-min', 'd-filling','WF', 'CN_M_sum', 'CN_O']
    obj = pd.DataFrame(db_dict, index = features_list ).T  # change the data into DataFrame format
    bar()

    obj.to_csv(f'{current_path}/NEW_BAND_CENTER_{sys.argv[1]}_ads_site_features_database_.csv')
    bar()

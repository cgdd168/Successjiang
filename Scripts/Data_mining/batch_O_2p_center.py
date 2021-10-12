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
import sys
import bisect
from ase.io import read, write
from scipy.integrate import simps
from alive_progress import alive_bar

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
def band_center(atom, orb, subfolder):
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

    bc_up = simps(y1 * x, x) / simps(y1, x)
    bc_down = simps(y2 * x, x) / simps(y2, x)
    bc_max = np.max([bc_up, bc_down])
    bc_min = np.min([bc_up, bc_down])
    dbc_avg = simps(y3 * x, x) / simps(y3, x)
    # band_width = np.sqrt( simps(y3 * x * x, x) / simps(y3, x) )
    print(subfolder, ' ', atom, f' atom {orb} band center [up, down, max, min, avg] is ', bc_up, bc_down, bc_max, bc_min, dbc_avg)
    return bc_up, bc_down, bc_max, bc_min, dbc_avg

# 执行计算任务的代码
global current_path
current_path = os.getcwd()
db_dict = {}
for foldername, subfolders, filenames in os.walk(current_path):
    with alive_bar(1 * len(subfolders), title='Total progress of the task') as bar:
        for i in range(len(subfolders)):
            cwd = os.chdir(subfolders[i])
            OUTPUT=[]                                     # 输出列表

            # 计算吸附位点O的p带中心细节
            p_up_center, p_down_center, p_max_center, p_min_center, p_band_center = band_center(atom=int(sys.argv[1]), orb='p', subfolder=subfolders[i])
            OUTPUT.append(p_up_center)
            OUTPUT.append(p_down_center)
            OUTPUT.append(p_max_center)
            OUTPUT.append(p_min_center)
            OUTPUT.append(p_band_center)
            bar()

            db_dict.update({subfolders[i]: OUTPUT})  # update the dict
            os.chdir(current_path)
    break

with alive_bar( 2, title='Writing to DATABASE') as bar:

    features_list = ['p_up_center', 'p_down_center', 'p_max_center', 'p_min_center', 'p_band_center']
    obj = pd.DataFrame(db_dict, index = features_list ).T  # change the data into DataFrame format
    bar()

    obj.to_csv(f'{current_path}/{sys.argv[1]}_ads_site_O_2p_details.csv')
    bar()

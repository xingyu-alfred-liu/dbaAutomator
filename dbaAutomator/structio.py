""" 
    Author: Xingyu (Alfred) Liu 
    Email: xingyu.alfred.liu@gmail.com
    Description:
        structio.py handles all io related functions
"""

import os
import sys
from pymatgen import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.xyz import XYZ
from ase.io import read
import numpy as np
import json

# singleMol is a dictionary, keys are atom index, values are pmg sites
# dataDir is the datapath, /dbaAutomator/data
def outputMolecule(singleMol, dataDir):
    molecule = Molecule([], [])
    singlemolpath = os.path.join(dataDir, 'singlemolecule')
    for siteIndex in singleMol.keys():
        molecule.append(str(singleMol[siteIndex].specie), singleMol[siteIndex].coords)
    xyzObj = XYZ(molecule)
    if "singleMol.xyz" in os.listdir(singlemolpath):
        decision = None
        print('You have one \'singleMol.xyz\' file inside the single molecule path.')
        print('This previous file will be overwritten.')
        while decision != 'Y' and decision != 'N':
            decision = input('Do you want to proceed? Y for yes, N for no.')
            if decision == 'Y':
                xyzObj.write_file(os.path.join(singlemolpath, 'singleMol.xyz'))
                print('The single molecule structure is saved under:', os.path.join(dataDir, 'singlemolecule.singleMol.xyz'))
            elif decision == 'N':
                print('The previous file is not changed. ')
                sys.exit()
            else:
                print('Not eligible response!!!\n')
    else:
        xyzObj.write_file(os.path.join(singlemolpath, 'singleMol.xyz'))
        print('The single molecule structure is saved under:', os.path.join(dataDir, 'singlemolecule/singleMol.xyz'))

# loadUnitCell handles both file and dir
# file has to be the DFT input unitcell file
# dir has to be the data dir
def loadUnitCell(path):
    if os.path.isfile(path):
        unitcell = read(path)
    elif os.path.isdir(path):
        unitcellpath = os.path.join(path, 'unitcell')
        filelist = os.listdir(unitcellpath)
        for name in filelist:
            try:
                unitcell = read(os.path.join(unitcellpath, name))
                if len(unitcell.get_positions()) != 0:
                    print('Please make sure this is the file you want to load:', os.path.join(unitcellpath, name))
                    break
            except:
                pass
    unitcell.set_pbc((True, True, True))
    pmgstruct = AseAtomsAdaptor.get_structure(unitcell)
    return pmgstruct

def loadSingleMol(path):
    singleMol = Molecule.from_file(os.path.join(path, 'singlemolecule/singleMol.xyz'))
    return singleMol

def loadChargeMatrix(struct, path):
    chargeMatrix = []
    try:
        # should also check if this charge file is correct
        with open(os.path.join(path, 'ACF.dat'), 'r') as infile:
            for values in infile:
                if len(values.split()) == 7 and values.split()[0] != '#':
                    chargeMatrix.append(values.split()[4])
                else:
                    pass
    except FileNotFoundError as fileerr:
        print(fileerr)
        print('Error !!!')
        print('Please make sure the Bader output is put into single Molecule Directory.')
        sys.exit()
    chargeMatrix = np.array(chargeMatrix)
    chargeMatrix = chargeMatrix.astype(float)
    if struct.num_sites != len(chargeMatrix):
        raise Exception("!!!Error!!!\nThe number of atoms does not match with molecule number")
    chargeSum = np.sum(chargeMatrix)
    chargeMatrix /= chargeSum
    # returned chargeMatrix provides the charge percentage and atom index
    return chargeMatrix

def outputHolePositions(holeSites, path):
    supercellPath = os.path.join(path, 'supercell')
    # json does not allow numpy,int64
    # convert values into list of float
    filename = 'holePositions.json'
    if filename in os.listdir(supercellPath):
        print('There is one json file from previous calculation')
        print('The old file will be rewritten')
        decision = 'A+'
        while decision != 'Y' and decision != 'N':
            decision = input('Do you want to proceed? \'Y\' for yes and \'N\' for no.')
            if decision == 'Y':
                pass
            elif decision == 'N':
                sys.exit('Exit now...')
            else:
                print('Please type in either Y or N.')
    tmpdict = dict()
    for key in holeSites.keys():
        tmpdict[int(key)] = list(holeSites[key])
    with open(supercellPath+'/holePositions.json', 'w') as file:
        file.write(json.dumps(tmpdict, indent=4))

def loadCubeCell(path):
    print('Loading supercell, this process might take minutes, please wait...')
    namelist = os.listdir(path)
    asestruct = None
    for name in namelist:
        if name.endswith('cube'):
            asestruct = read(os.path.join(path, name))
            break
    if asestruct == None:
        print('Error!!!')
        print('There is no cube file under /data/supercell.')
        print('Please make sure the exciton wavefunction calculation output\
            is put under the folder \'supercell\'.')
        sys.exit()
    else:
        pmgstruct = AseAtomsAdaptor.get_structure(asestruct)
        return pmgstruct

def loadHolePositions(path):
    supercellpath = os.path.join(path, "supercell")
    try:
        with open(os.path.join(supercellpath, "holePositions.json"), 'r') as jsonin:
            data = json.load(jsonin)
            return data
    except FileNotFoundError as fileerr:
        print(fileerr)
        sys.exit("Please make sure file holePositions.json is put under /data/supercell.")

def printChargeShare(a, b, c, convThreshold):
    print('The charge share for direction a is:', "{:0.2f}".format(a*100), "%.")
    print('The charge share for direction b is:', "{:0.2f}".format(b*100), "%.")
    print('The charge share for direction c is:', "{:0.2f}".format(c*100), "%.")
    if a > convThreshold:
        print('!!! Warning !!!')
        print('Boundary atoms in direction a occupies more than', "{:0.2f}".format(convThreshold*100), "%", "charge.")
    if b > convThreshold:
        print('!!! Warning !!!')
        print('Boundary atoms in direction b occupies more than', "{:0.2f}".format(convThreshold*100), "%", "charge.")
    if c > convThreshold:
        print('!!! Warning !!!')
        print('Boundary atoms in direction c occupies more than', "{:0.2f}".format(convThreshold*100), "%", "charge.")

# this function will create a list of directories
# name them with the charge site index of the single molecule
# and place the plot_xct input files there
def createPlotxctInput(path, holeSites, fineGrid):
    dbapath = os.path.join(path, 'dba')
    for key in holeSites.keys():
        plotxctpath = os.path.join(dbapath, str(key))
        os.system('mkdir '+plotxctpath)
        with open(os.path.join(plotxctpath, 'plotxct.inp'), 'w') as outfile:
            outfile.write("# Index of state to be plotted, as it appears in eigenvectors\nplot_state 1\n")
            outfile.write("# Size of supercell\n")
            outfile.write("supercell_size  ")
            for grid in fineGrid:
                outfile.write(str(grid)+'  ')
            outfile.write('\n')
            outfile.write("# coordinates of hole, in units of lattice vectors\n")
            outfile.write('hole_position   ')
            for value in holeSites[key]:
                outfile.write(str(value)+'  ')

def loadPlotxct(path):
    try:
        with open(os.path.join(path, "plotxct.inp"), 'r') as f:
            for line in f:
                if "hole_position" in line:
                    tmpline = line.split()
                    holePosition = tmpline[1:]
        for i, val in enumerate(holePosition):
            holePosition[i] = float(val)
    except FileNotFoundError as fileerr:
        print(fileerr)
        print('Please make sure your plotxct.inp is put under:', path)
        sys.exit()
    return holePosition

# write the result of dba into a folder
def writedbaResult(path, chargeshare, chargetransfer):
    supercellPath = os.path.join(path, 'supercell')
    # json does not allow numpy,int64
    # convert values into list of float
    filename = 'dba.txt'
    if filename in os.listdir(supercellPath):
        print('There is one dba.txt from previous calculation')
        print('The old file will be rewritten')
        decision = None
        while decision != 'Y' and decision != 'N':
            decision = input('Do you want to proceed? \'Y\' for yes and \'N\' for no.')
            if decision == 'Y':
                pass
            elif decision == 'N':
                sys.exit('Exit now...')
            else:
                print('Please type in either Y or N.\n')
    with open(os.path.join(supercellPath, filename), 'w') as file:
        for hole in chargeshare.keys():
            file.write("The charge transfer for hole "+hole+" is: "+"{:0.2f}".format((1-chargeshare[hole])*100)+" %.\n")
        file.write("The total charge transfer is: "+"{:0.2f}".format((1-chargetransfer)*100)+" %.\n")
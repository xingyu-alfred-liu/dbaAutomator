import os
import sys
from pymatgen import Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.xyz import XYZ
from ase.io import read
import numpy as np
import json

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
                print('The single molecule structure and corresponding index in the supercell is saved under \'/data/singlemolecule\'')
                print('It\'s named as \'singleMol.xyz\'.')
            elif decision == 'N':
                print('The previous file is not changed. ')
                sys.exit()
            else:
                print('Not eligible response!!!\n')
    else:
        xyzObj.write_file(os.path.join(singlemolpath, 'singleMol.xyz'))
        print('The single molecule structure and corresponding index in the supercell is saved under \'/data/singlemolecule\'')
        print('It\'s named as \'singleMol.xyz\'.')

def loadUnitCell(path):
    unitcellpath = path+'/unitcell'
    fileNum = len(os.listdir(unitcellpath))
    if fileNum == 0:
        print('Error!!!')
        print('Please include unit cell information under /data/unitcell')
        return 0
    else:
        filelist = os.listdir(unitcellpath)
        if 'kgrid.in' in filelist:
            filelist.remove('kgrid.in')
        for filename in filelist:
            if filename.startswith('.'):
                filelist.remove(filename)
        unitcell = read(os.path.join(unitcellpath, filelist[0]))
        if unitcell == 0:
            print('There is no readable input file in /data/unitcell')
            print('Please check the instructions or change your unit cell format')
            return 0
        else:
            unitcell.set_pbc((True, True, True))
            pmgobj = AseAtomsAdaptor()
            pmgstruct = pmgobj.get_structure(unitcell)
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
                if len(values.split()) == 7 and values.split()[0].isdigit():
                    chargeMatrix.append(values.split())
                else:
                    pass
    except FileNotFoundError as fileerr:
        print(fileerr)
        print('!!! Error !!!')
        print('Please make sure the Bader output is put into single Molecule Directory.')
        return 0
    chargeMatrix = np.array(chargeMatrix)
    chargeMatrix = chargeMatrix.astype(float)
    if struct.num_sites != len(chargeMatrix):
        raise Exception("!!!Error!!!\nThe number of atoms does not match with molecule number")
    chargeSum = np.sum(chargeMatrix, axis=0)
    chargeMatrix[:, 4] /= chargeSum[4]
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
            decision = input('Do you want to proceed? Y for \'yes\' and N for \'no\'.\n')
            if decision == 'Y':
                pass
            elif decision == 'N':
                return 0
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
        pmgobj = AseAtomsAdaptor()
        pmgstruct = pmgobj.get_structure(asestruct)
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
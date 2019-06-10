"""
This script gives the hole positions w.r.t supercell generated before
To locate the hole position, we need the single molecule HOMO charge,
bader analysis result (ACF.dat) and the unit cell. 

"""
from filepath import *
from ase.io import read, write
import os
import numpy as np

def loadSingleMol(singleMolPath):
    singleMol = read(singleMolPath+'/singleMol.xyz')
    return singleMol

def getChargeMatrix(singleMol, singleMolPath):
    if not 'ACF.dat' in os.listdir(singleMolPath):
        print('!!! Error !!!')
        print('Please make sure the Bader output is put into single Molecule Directory.')
        return 0
    chargeMatrix = []
    # should also check if this charge file is correct
    with open(singleMolPath+'/ACF.dat', 'r') as infile:
        for values in infile:
            if len(values.split()) == 7 and values.split()[0].isdigit():
                chargeMatrix.append(values.split())
            else:
                pass
    chargeMatrix = np.array(chargeMatrix)
    chargeMatrix = chargeMatrix.astype(float)
    if singleMol.get_number_of_atoms() != len(chargeMatrix):
        print('!!! Error !!!')
        print('The number of atoms does not match with molecule number')
        return 0
    chargeSum = np.sum(chargeMatrix, axis=0)
    chargeMatrix[:, 4] /= chargeSum[4]
    # returned chargeMatrix provides the charge percentage and atom index
    return chargeMatrix

if __name__ == "__main__":
    singleMol = loadSingleMol(singleMolPath)
    chargeMatrix = getChargeMatrix(singleMol, singleMolPath)
    
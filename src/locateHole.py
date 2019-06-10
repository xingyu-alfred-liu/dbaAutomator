"""
This script gives the hole positions w.r.t supercell generated before
To locate the hole position, we need the single molecule HOMO charge,
bader analysis result (ACF.dat) and the unit cell. 

"""
from filepath import *
import os
import numpy as np
from selectSingleMol import constructSuperCell, getBondDict
from pymatgen import Molecule
from bondLengthRef import bondCutoff

def loadSingleMol(singleMolPath):
    singleMol = Molecule.from_file(singleMolPath+'/singleMol.xyz')
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
    if singleMol.num_sites != len(chargeMatrix):
        print('!!! Error !!!')
        print('The number of atoms does not match with molecule number')
        return 0
    chargeSum = np.sum(chargeMatrix, axis=0)
    chargeMatrix[:, 4] /= chargeSum[4]
    # returned chargeMatrix provides the charge percentage and atom index
    return chargeMatrix

# !!! important !!!
# this function sets the charge percentage threshold as 1%
def getHolePositions(chargeMatrix, singleMol, supercell, bondDict, chargeThreshold=0.1):
    bondlength = 0
    for key in bondDict:
        if bondDict[key] >= bondlength:
            bondlength = bondDict[key]
    print(supercell.sites[0].frac_coords)
    chargeIndex = np.where(chargeMatrix[:, 4] > chargeThreshold)[0]
    holePositions = []
    for charindex in chargeIndex:
        chargeSite = singleMol.sites[charindex]
        neighborSites = singleMol.get_neighbors(singleMol.sites[charindex], bondlength)
    return None

if __name__ == "__main__":
    singleMol = loadSingleMol(singleMolPath)
    chargeMatrix = getChargeMatrix(singleMol, singleMolPath)
    supercell = constructSuperCell(QEinputPath, fineGrid, fineGridpath)
    # find out the fractional coordinates of hole positions
    bondDict = getBondDict(supercell, bondCutoff)
    holeSites = getHolePositions(chargeMatrix, singleMol, supercell, bondDict, chargeThreshold)
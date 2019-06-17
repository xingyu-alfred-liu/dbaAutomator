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
from copy import deepcopy
import math
from itertools import filterfalse
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor

def loadSingleMol(singleMolPath):
    singleMol = Molecule.from_file(singleMolPath+'/singleMol.xyz')
    return singleMol

def getChargeMatrix(struct, path):
    if not 'ACF.dat' in os.listdir(path):
        print('!!! Error !!!')
        print('Please make sure the Bader output is put into single Molecule Directory.')
        return 0
    chargeMatrix = []
    # should also check if this charge file is correct
    with open(path+'/ACF.dat', 'r') as infile:
        for values in infile:
            if len(values.split()) == 7 and values.split()[0].isdigit():
                chargeMatrix.append(values.split())
            else:
                pass
    chargeMatrix = np.array(chargeMatrix)
    chargeMatrix = chargeMatrix.astype(float)
    if struct.num_sites != len(chargeMatrix):
        print('!!! Error !!!')
        print('The number of atoms does not match with molecule number')
        return 0
    chargeSum = np.sum(chargeMatrix, axis=0)
    chargeMatrix[:, 4] /= chargeSum[4]
    # returned chargeMatrix provides the charge percentage and atom index
    return chargeMatrix

def calNormalVector(p1, p2, p3):
    # generate the vector
    vector = [0, 0, 0]
    vector[0] = (p2[1]-p1[1])*(p3[2]-p1[2]) - (p2[2]-p1[2])*(p3[1]-p1[1])
    vector[1] = (p2[2]-p1[2])*(p3[0]-p1[0]) - (p2[0]-p1[0])*(p3[2]-p1[2])
    vector[2] = (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0])
    # normalize it
    sigma = vector[0]**2 + vector[1]**2 + vector[2]**2
    sigma = math.sqrt(sigma)
    for i in range(3):
        vector[i] = vector[i] / sigma
    return vector

def findHole(unitcell, twoNeighbors, chargeSite):
    point1 = deepcopy(twoNeighbors[0][0].coords)
    point2 = deepcopy(twoNeighbors[1][0].coords)
    point3 = deepcopy(chargeSite.coords)
    normalVec = calNormalVector(point1, point2, point3)
    shift = -0.8
    holePosition = [0, 0, 0]
    for i in range(3):
        holePosition[i] = chargeSite.coords[i] + normalVec[i]*shift
    unitcell.append('He', holePosition, coords_are_cartesian=True)
    print('charge site', chargeSite)
    print('hole position', holePosition)
    return unitcell.sites[-1].frac_coords

# !!! important !!!
# this function sets the charge percentage threshold as 1%
def getHolePositions(chargeMatrix, singleMol, unitcell, bondDict, chargeThreshold=0.01):
    bondlength = 0
    for key in bondDict:
        if bondDict[key] >= bondlength:
            bondlength = bondDict[key]
    chargeIndex = np.where(chargeMatrix[:, 4] > chargeThreshold)[0]
    holePositions = dict()
    for charindex in chargeIndex:
        chargeSite = singleMol.sites[charindex]
        neighborSites = singleMol.get_neighbors(chargeSite, bondlength)
        twoNeighbors = []
        # delete the neighbors which are H
        print()
        neighborSites[:] = filterfalse(lambda x: str(x[0].specie) == 'H', neighborSites)
        # print('charge site')
        # print(chargeSite)
        # print('neighbor site')
        # print(neighborSites)
        # control the length of the twoNeighbors list smaller or equal than 2
        # delete every
        for neighbor in neighborSites:
            if len(twoNeighbors) < 2:
                twoNeighbors += [neighbor]
            else:
                twoNeighbors += [neighbor]
                tmpBondLength = 0
                for site in twoNeighbors:
                    if site[1] > tmpBondLength:
                        tmpBondLength = site[1]
                        removeSite = site
                twoNeighbors.remove(removeSite)
        holePositions[charindex] = findHole(unitcell, twoNeighbors, chargeSite)
    print()
    print('hole positions:')
    for key in holePositions.keys():
        print(key, ':', holePositions[key])
    return holePositions

def getUnitCell(QEinputPath):
    asecell = read(QEinputPath)
    mgobj = AseAtomsAdaptor()
    mgcell = mgobj.get_structure(asecell)
    return mgcell

if __name__ == "__main__":
    singleMol = loadSingleMol(singleMolPath)
    chargeMatrix = getChargeMatrix(singleMol, singleMolPath)
    #supercell = constructSuperCell(QEinputPath, fineGrid, fineGridpath)
    unitcell = getUnitCell(QEinputPath)
    # find out the fractional coordinates of hole positions
    bondDict = getBondDict(unitcell, bondCutoff)
    #holeSites = getHolePositions(chargeMatrix, singleMol, supercell, bondDict, chargeThreshold)
    holeSites = getHolePositions(chargeMatrix, singleMol, unitcell, bondDict, chargeThreshold)
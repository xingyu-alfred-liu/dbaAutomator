from filepath import *
from pymatgen import Molecule
from locateHole import loadSingleMol, getChargeMatrix
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
import os
import numpy as np

def getSuperCell(path):
    namelist = os.listdir(path)
    for name in namelist:
        if name.endswith('cube'):
            cubefile = name
            break
    asestruct = read(path+'/'+cubefile)
    mgobj = AseAtomsAdaptor()
    mgstruct = mgobj.get_structure(asestruct)
    return mgstruct

# find out the index for each site in the supercell
# return the dict, key is singlemolecule index
# value is the supercell index
def getMoleculeIndex(singleMol, cubecell, threshold = 0.01):
    molIndex = dict()
    for i, molsite in enumerate(singleMol.sites):
        for j, cellsite in enumerate(cubecell.sites):
            if (np.linalg.norm(molsite.coords-cellsite.coords)<threshold) and (str(cellsite.specie)==str(molsite.specie)):
                molIndex[i] = j
            continue
    print('The length of the hole positioned molecule:', len(molIndex.keys()))
    return molIndex

def getMolShare(chargeMatrix, molIndex):
    molshare = 0
    for key in molIndex.keys():
        molshare += chargeMatrix[molIndex[key]][4]
    return molshare

if __name__ == "__main__":
    singleMol = loadSingleMol(singleMolPath)
    cubeSuperCell = getSuperCell(supercellPath)
    molIndex = getMoleculeIndex(singleMol, cubeSuperCell)
    chargeMatrix = getChargeMatrix(cubeSuperCell, supercellPath)
    molShare = getMolShare(chargeMatrix, molIndex)
    print('The charge share for the hole positioned molecule is:', "{:0.2f}".format(molShare*100), "%")
from ase.io import read, write
from ase.build import make_supercell
from filepath import *
import linecache
from copy import deepcopy
import numpy as np
# ase has some problem in its make supercell function so we have
# switch to pymatgen
from pymatgen.io.ase import AseAtomsAdaptor
from bondLengthRef import bondCutoff
from pymatgen import Molecule

# getBondDict give the bond cutoff dict for supercell
def getBondDict(supercell, bondCutoff):
    bondDict = dict()
    speciesList = list(set(supercell.species))
    for i in range(len(speciesList)):
        speciesList[i] = str(speciesList[i])
    for ele in speciesList:
        for pairs in bondCutoff.keys():
            if (ele in pairs) and (pairs not in bondDict):
                bondDict[pairs] = bondCutoff[pairs]
    duplicate = dict()
    for pairs in bondDict.keys():
        (a, b) = pairs
        if (b, a) not in bondDict:
            duplicate[(b, a)] = bondDict[pairs]
    bondDict.update(duplicate)
    return bondDict

def getFineGrid(path):
    filein = linecache.getlines(path)
    fineGrid = filein[0].split()
    for i in range(len(fineGrid)):
        fineGrid[i] = int(fineGrid[i])
    return fineGrid

def constructSuperCell(QEinputPath, fineGrid):
    # ase has some problem in its super cell constructing function
    # here switch to pymatgen
    structure = read(QEinputPath)
    structure.set_pbc((True, True, True))
    pymatgenObject = AseAtomsAdaptor()
    pymatgenStruct = pymatgenObject.get_structure(structure)
    pymatgenStruct.make_supercell(fineGrid)
    # pymatgenStruct.to(filename='pymatgenSuperCell.cif')
    return pymatgenStruct

def getSingleMol(supercell, bondDict):
    # find site in the middle
    dist = 1
    middle = [0.5, 0.5, 0.5]
    for site in supercell.sites:
        if np.linalg.norm(middle-site.frac_coords) < dist:
            dist = np.linalg.norm(middle-site.frac_coords)
            middleSite = site
    print(site.frac_coords)
    print(site.coords)
    print(site.specie)
    site1 = supercell.sites[0]
    print(supercell.get_neighbors(site1, 1))
    print('The site closest to the middle is', middleSite)
    MolEleList = []
    MolEleCoords = []
    return None

if __name__ == "__main__":

    print('To users: please make sure you type in correct path for necessary files in \"filepath.py\"')
    
    # decide the fine grid
    if fineGrid == []:
        fineGrid = getFineGrid(fineGridpath)
    # read in the structure file and construct the super cell
    supercell = constructSuperCell(QEinputPath, fineGrid)
    # get the single molecule from the super cell in the middle
    bondDict = getBondDict(supercell, bondCutoff)
    singleMol = getSingleMol(supercell, bondDict)
from ase.io import read, write
from ase.build import make_supercell
from filepath import *
import linecache
from copy import deepcopy
import numpy as np
# ase has some problem in its make supercell function so we have
# switch to pymatgen
from pymatgen.io.ase import AseAtomsAdaptor

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

def getSingleMol(supercell):
    # find site in the middle
    dist = 1
    middle = [0.5, 0.5, 0.5]
    for site in supercell.sites:
        if np.linalg.norm(middle-site.frac_coords) < dist:
            dist = np.linalg.norm(middle-site.frac_coords)
            middleSite = site
    print('The site closest to the middle is', middleSite)
    return None

if __name__ == "__main__":

    print('To users: please make sure you type in correct path for necessary files in \"filepath.py\"')
    
    # decide the fine grid
    if fineGrid == []:
        fineGrid = getFineGrid(fineGridpath)
    # read in the structure file and construct the super cell
    supercell = constructSuperCell(QEinputPath, fineGrid)
    # get the single molecule from the super cell in the middle
    singleMol = getSingleMol(supercell)
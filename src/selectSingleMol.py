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
from pymatgen.io.xyz import XYZ
import os

# getBondDict give the bond cutoff dict for supercell
def getBondDict(supercell, bondCutoff):
    bondDict = dict()
    speciesList = list(set(supercell.species))
    for i in range(len(speciesList)):
        speciesList[i] = str(speciesList[i])
    for a in speciesList:
        for b in speciesList:
            if (a, b) in bondCutoff.keys():
                bondDict[(a, b)] = bondCutoff[(a, b)]
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

def constructSuperCell(QEinputPath, fineGrid, fineGridpath):
    # decide the fine grid
    if fineGrid == []:
        fineGrid = getFineGrid(fineGridpath)
    # ase has some problem in its super cell constructing function
    # here switch to pymatgen
    structure = read(QEinputPath)
    structure.set_pbc((True, True, True))
    pymatgenObject = AseAtomsAdaptor()
    pymatgenStruct = pymatgenObject.get_structure(structure)
    pymatgenStruct.make_supercell(fineGrid)
    # pymatgenStruct.to(filename='pymatgenSuperCell.cif')
    return pymatgenStruct

def getSingleMol(supercell, middleSite, bondDict, middleSiteIndex):
    candidates = [middleSite]
    candidatesIndex = [middleSiteIndex]
    tmpSites = []
    singleMol = dict()
    while candidates != []:
        for site in candidates:
            bondCutoffMax = 0
            for pair in bondDict.keys():
                if (bondDict[pair] >= bondCutoffMax) and (str(site.specie) in pair):
                    bondCutoffMax = bondDict[pair]
            allNeighbors = supercell.get_neighbors(site, bondCutoffMax, include_index=True)
            for neighbor in allNeighbors:
                sitePair = (str(site.specie), str(neighbor[0].specie))
                if neighbor[1] <= bondDict[sitePair]:
                    tmpSites += [neighbor]
        tmpSites = list(set(tmpSites))
        for i, site in enumerate(candidates):
            singleMol[candidatesIndex[i]] = site
        candidates = []
        candidatesIndex = []
        for site in tmpSites:
            # check if the site index is in the single molecule index
            if site[2] not in singleMol.keys():
                candidates += [site[0]]
                candidatesIndex += [site[2]]
        print('Number of atoms in this molecule:')
        print(len(singleMol))
        tmpSites = []
    return singleMol

def getCentralSingleMol(supercell, bondDict):
    # find site in the middle
    dist = 1
    middle = [0.5, 0.5, 0.5]
    for i, site in enumerate(supercell.sites):
        if np.linalg.norm(middle-site.frac_coords) < dist:
            dist = np.linalg.norm(middle-site.frac_coords)
            middleSite = site
            middleSiteIndex = i
    # print(site.frac_coords)
    # print(site.coords)
    # print(site.specie)
    print('The site closest to the middle is', middleSite)
    print('The corresponding site index is', middleSiteIndex)
    # pick up all the atom pairs within bond van der waals distance
    # centralSingleMol is the list of sites which belong to the central molecule
    centralSingleMol = getSingleMol(supercell, middleSite, bondDict, middleSiteIndex)
    return centralSingleMol

def outputMolecule(singleMol, dataDir):
    molecule = Molecule([], [])
    moleculeIndex = []
    for siteIndex in singleMol.keys():
        moleculeIndex.append(siteIndex)
        molecule.append(str(singleMol[siteIndex].specie), singleMol[siteIndex].coords)
    moleculeIndex = np.array(moleculeIndex)
    moleculeIndex = moleculeIndex.astype(int)
    xyzObj = XYZ(molecule)
    os.chdir(dataDir)
    os.system('mkdir singleMolecule')
    xyzObj.write_file(dataDir+'/singleMol.xyz')
    np.savetxt(dataDir+'/singleMolIndex.txt', moleculeIndex, delimiter=',', fmt="%d")
    print('The single molecule structure and corresponding index in the supercell is saved under \'singleMol\'')

if __name__ == "__main__":

    print('To users: please make sure you type in correct path for necessary files in \"filepath.py\"')
    
    # read in the structure file and construct the super cell
    supercell = constructSuperCell(QEinputPath, fineGrid, fineGridpath)
    # get the single molecule from the super cell in the middle
    bondDict = getBondDict(supercell, bondCutoff)
    singleMol = getCentralSingleMol(supercell, bondDict)
    # otuput the single Molecule file
    outputMolecule(singleMol, singleMolPath)
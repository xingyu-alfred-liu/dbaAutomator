""" 
    utility.py provide the io functions, including:

    - loadUnitCell(path)
    path is the dataPath, dataPath needs to be defined
    in the files.py
    
"""

from ase.io import read, write
import os
from pymatgen.io.ase import AseAtomsAdaptor
import linecache
import numpy as np
from copy import deepcopy
from pymatgen import Molecule
from pymatgen.io.xyz import XYZ

def outputMolecule(singleMol, dataDir):
    molecule = Molecule([], [])
    #moleculeIndex = []
    for siteIndex in singleMol.keys():
        #moleculeIndex.append(siteIndex)
        molecule.append(str(singleMol[siteIndex].specie), singleMol[siteIndex].coords)
    #moleculeIndex = np.array(moleculeIndex)
    #moleculeIndex = moleculeIndex.astype(int)
    xyzObj = XYZ(molecule)
    xyzObj.write_file(dataDir+'/singlemolecule/singleMol.xyz')
    #np.savetxt(dataDir+'/singleMolIndex.txt', moleculeIndex, delimiter=',', fmt="%d")
    print('The single molecule structure and corresponding index in the supercell is saved under \'/data/singlemolecule\'')

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
        # ugly solution
        tmpSiteslist = list()
        for tmpsite in tmpSites:
            if tmpsite not in tmpSiteslist:
                tmpSiteslist += [tmpsite]
        tmpSites = deepcopy(tmpSiteslist)
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
    print('The site closest to the middle is', middleSite)
    print('The corresponding site index is', middleSiteIndex)
    # pick up all the atom pairs within bond van der waals distance
    # centralSingleMol is the list of sites which belong to the central molecule
    centralSingleMol = getSingleMol(supercell, middleSite, bondDict, middleSiteIndex)
    return centralSingleMol

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
        unitcell = read(unitcellpath+'/'+filelist[0])
        if unitcell == 0:
            print('There is no readable input file in /data/unitcell')
            print('Please check the instructions or change your unit cell format')
            return 0
        else:
            unitcell.set_pbc((True, True, True))
            pmgobj = AseAtomsAdaptor()
            pmgstruct = pmgobj.get_structure(unitcell)
    return pmgstruct

def getSuperCell(path, unitcell, finegrid):
    unitcellpath = path+'/unitcell'
    if finegrid != []:
        unitcell.make_supercell(finegrid)
    else:
        print('No definitions found in fineGrid, now trying kgrid.in')
        if 'kgrid.in' in os.listdir(unitcellpath):
            finegrid = getFineGrid(unitcellpath+'/kgrid.in')
            print('Fine grid found in kgrid.in is:', finegrid)
        else:
            print('Error!!!')
            print('There is no definition of fineGrid and no \'kgrid.in\' found in unitcell folder, please check your settings.')
            return 0
        unitcell.make_supercell(finegrid)
    return unitcell
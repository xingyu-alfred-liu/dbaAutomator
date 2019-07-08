""" 
    utility.py provide all functions:

    - loadUnitCell(path)
    path is the dataPath, dataPath needs to be defined
    in the files.py
    
"""

import os
import sys
import numpy as np
from copy import deepcopy
from pymatgen import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from itertools import filterfalse
import math

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

def getSuperCell(unitcell, finegrid):
    if finegrid != []:
        print()
        print('Reminder: Please make sure you type the correct fine grid')
        print()
        unitcell.make_supercell(finegrid)
    else:
        sys.exit('No definitions found in finegrid, please make sure you input it when define automator.\n')
    return unitcell

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
        neighborSites[:] = filterfalse(lambda x: str(x[0].specie) == 'H', neighborSites)
        # control the length of the twoNeighbors list smaller or equal than 2
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
    return holePositions

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
    return unitcell.sites[-1].frac_coords

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

# this function will create a list of directories
# name them with the charge site index of the single molecule
# and place the plot_xct input files there
def createPlotxctInput(path, holeSites, fineGrid):
    dbapath = os.path.join(path, 'dba')
    for key in holeSites.keys():
        plotxctpath = os.path.join(dbapath, str(key))
        os.system('mkdir '+plotxctpath)
        with open(os.path.join(plotxctpath, 'plotxct.inp'), 'w') as outfile:
            outfile.write("Index of state to be plotted, as it appears in eigenvectors\n\
                           plot_state 1\n")
            outfile.write("# Size of supercell\n")
            outfile.write("supercell_size  ")
            for grid in fineGrid:
                outfile.write(str(grid)+'  ')
            outfile.write("# coordinates of hole, in units of lattice vectors\n")
            outfile.write('hole_position   ')
            for value in holeSites[key]:
                outfile.write(str(value)+'  ')

def getMoleculeIndex(singleMol, cubecell, threshold = 0.01):
    molIndex = dict()
    for i, molsite in enumerate(singleMol.sites):
        for j, cellsite in enumerate(cubecell.sites):
            if (np.linalg.norm(molsite.coords-cellsite.coords)<threshold) and (str(cellsite.specie)==str(molsite.specie)):
                molIndex[i] = j
            continue
    print('The length of the hole positioned molecule:', len(molIndex.keys()))
    if len(molIndex.keys()) != singleMol.num_sites:
        print('Error!!!')
        print('The number of atoms within a single molecule found in the cube file\
            is not the same as the one stored under /data/singlemolecule')
        print('Please make sure there is no change of the output single \
            molecule from step 1')
        return 0
    else:
        return molIndex

def getMolShare(chargeMatrix, molIndex):
    molshare = 0
    for key in molIndex.keys():
        molshare += chargeMatrix[molIndex[key]][4]
    return molshare
    
def getXctPath(path, checklist):
    filelist = os.listdir(path)
    if "ACF.dat" not in filelist:
        for name in filelist:
            if os.path.isdir(os.path.join(path, name)):
                checklist = getXctPath(os.path.join(path, name), checklist)
    else:
        if not any(name.endswith("cube") for name in filelist):
            print()
            print("Warning!!!")
            print('In folder', path, 'we found ACF.dat but no .cube file.')
            print('Please check if you put all necessary output there.')
        else:
            checklist += [path]
    return checklist

def getPrimitiveCell(supercell):
    pmganalyzer = SpacegroupAnalyzer(supercell)
    primitive = pmganalyzer.find_primitive()
    return primitive

# getAllMols returns a list of pmg Molecule objs, including all fragments inside the supercell
def getAllMols(supercell, bondDict):
    # molList is the list of pymatgen Molecule object
    # tmpMol is pyatgen Molecule object
    molList = []
    while len(supercell.sites) != 0:
        # function getCentralSingleMol returns a dictionary
        # the keys are the index of each sites in supercell
        print('length of supercell:', len(supercell.sites))
        molsites = Molecule([], [])
        molindex = list()
        tmpMol = getCentralSingleMol(supercell, bondDict)
        for siteIndex in tmpMol.keys():
            molsites.append(str(tmpMol[siteIndex].specie), tmpMol[siteIndex].coords)
            molindex.append(siteIndex)
        supercell.remove_sites(molindex)
        molList += [molsites]
    return molList

# take a list of Molecule objects and return the smallest intermolecular distances
def getInterMolLen(molslist):
    mollen = 0
    for mol in molslist:
        if len(mol.sites) > mollen:
            mollen = len(mol.sites)
    for mol in molslist[:]:
        if len(mol.sites) < mollen:
            molslist.remove(mol)
    if len(molslist) >= 2:
        # choose two random molecules and get an initial intermolecular distance
        interLen = np.linalg.norm(molslist[0].center_of_mass - molslist[1].center_of_mass)
        for i, mola in enumerate(molslist):
            for j, molb in enumerate(molslist):
                # make sure not counting the difference between the same molecules
                if (i != j) and (np.linalg.norm(mola.center_of_mass-molb.center_of_mass) < interLen):
                    interLen = np.linalg.norm(mola.center_of_mass-molb.center_of_mass)
        return interLen
    else:
        raise Exception('Error!!! There are less than two complete molecules inside constructed supercell.\n')

# take a pmg structure and convergence length
# return the atom index in three dimensions
def getAtomIndex(struct, convrange):
    dirabc = list()
    for d in range(3):
        cutoff = convrange / struct.lattice.abc[d]
        dirtmp = np.where(np.logical_or(struct.frac_coords[:, d] < cutoff, struct.frac_coords[:, d] > (1-cutoff)))[0]
        dirabc += [dirtmp]
    return dirabc[0], dirabc[1], dirabc[2]

# take the indices and chargematrix, return
# the charge share for indicated indices
def getChargeShare(indices, chargematrix):
    return np.sum(chargematrix[indices], axis=0)[4]
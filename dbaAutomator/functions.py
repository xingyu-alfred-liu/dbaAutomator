""" 
    Author: Xingyu (Alfred) Liu 
    Email: xingyu.alfred.liu@gmail.com
    Description:
        functions.py provide all functions, input 
        and feature explanations can be found 
        around the functions.
"""

import os
import sys
import numpy as np
from copy import deepcopy
from pymatgen import Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.xyz import XYZ
from itertools import filterfalse
from shutil import copy2
import math

# getSingleMol finds a singlemolecule with provided crystal and the
# starting pymatgen site and index
# return: dict, key is index of the site, value is the site
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
        print('Number of atoms found in this molecule:', len(singleMol))
        tmpSites = []
    return singleMol

# getCentralSingleMol finds the atom closest to the middle of the crystal
# and use getSingleMol return the dictionay using index as key and sites 
# as values
# return: same as getSingleMol
def getCentralSingleMol(supercell, bondDict, middle=[0.5, 0.5, 0.5]):
    # check if the frac_coord is smaller than zero
    # if smaller than zero, change the middle site coords to negative
    if min(supercell.frac_coords[:, 0]) < 0:
        for i, val in enumerate(middle):
            middle[i] = val * (-1)
    # find site in the middle
    dist = 1
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

# getBondDict use bondCutoff reference and input structure
# to decide which kind of bonding should be selected
# return: dictionary, using atom pair as key, and the vdW distance as value
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

# getSuperCell makes a super using pymatgen make_supercell
# return: pymatgen Structure object
def getSuperCell(unitcell, finegrid):
    if finegrid != []:
        print('Reminder: Please make sure you type in the correct fine grid.')
        unitcell.make_supercell(finegrid)
    else:
        sys.exit('No definitions found in finegrid, please make sure you input it when define automator.\n')
    return unitcell

# !!! important !!!
# this function sets the charge percentage threshold as 1%
# getHolePositions finds out the hole positions with given bader output (chargeMatrix) and 
# structure
# return: dictionary, with charge site index as key and hole position as value
def getHolePositions(chargeMatrix, singleMol, unitcell, bondDict, chargeThreshold, holeAtomDist):
    bondlength = 0
    for key in bondDict:
        if bondDict[key] >= bondlength:
            bondlength = bondDict[key]
    chargeIndex = np.where(chargeMatrix > chargeThreshold)[0]
    holeProb = 0
    for charindex in chargeIndex:
        if chargeMatrix[charindex] > holeProb:
            holeHighProbIndex = charindex
            holeProb = chargeMatrix[charindex]
    print('The hole with index: %s has the highest probability.\n' %str(holeHighProbIndex))
    print('Please consider converge the exciton wave function calculated with this hole first.\n')
    print('You can use the checker object for convergence check.\n')
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
        holePositions[charindex] = findHole(unitcell, twoNeighbors, chargeSite, holeAtomDist)
    return holePositions

# find out hole position with given structure and charge site and corresponding two
# neighbors
# return hole positions fractional coords, using unitcell lattice vectors
def findHole(unitcell, twoNeighbors, chargeSite, holeAtomDist):
    point1 = deepcopy(twoNeighbors[0][0].coords)
    point2 = deepcopy(twoNeighbors[1][0].coords)
    point3 = deepcopy(chargeSite.coords)
    normalVec = calNormalVector(point1, point2, point3)
    shift = holeAtomDist
    holePosition = [0, 0, 0]
    for i in range(3):
        holePosition[i] = chargeSite.coords[i] + normalVec[i]*shift
    tmpcell = unitcell.copy()
    tmpcell.append('He', holePosition, coords_are_cartesian=True)
    return tmpcell.sites[-1].frac_coords

# calculate the normal vector with given three coords
# return: list, indicating the normal vector
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

# finds the site index with given search range from a structure
# return: dict, molecular site index as key and structure site index as value
def getMoleculeIndex(singleMol, cubecell, threshold = 0.01):
    molIndex = dict()
    for i, molsite in enumerate(singleMol.sites):
        for j, cellsite in enumerate(cubecell.sites):
            if (np.linalg.norm(molsite.coords-cellsite.coords)<threshold) and (str(cellsite.specie)==str(molsite.specie)):
                molIndex[i] = j
            continue
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
        molshare += chargeMatrix[molIndex[key]]
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
    return np.sum(chargematrix[indices])

def checkDataFolder(path):
    folderlist = ['dba', 'singlemolecule', 'supercell', 'unitcell']
    for foldername in folderlist:
        try:
            os.mkdir(os.path.join(path, foldername))
            print('Directory', foldername, 'created.')
        except FileExistsError:
            print('Directory', foldername, 'already exists.')

def copyInput(file, path):
    unitcellpath = os.path.join(path, 'unitcell')
    print('Copy unitcell file to', unitcellpath, '.')
    copy2(file, unitcellpath)

def getMPC(supercell, finegrid, molslist):
    # get the number of unitcells used to build the supercell
    # get the number of atoms per cell
    supercellsize = finegrid[0] * finegrid[1] * finegrid[2]
    atomPerCell = len(supercell.sites) / supercellsize
    if len(molslist) == 0:
        sys.exit('There is no complete molecules founded in molslist...')
    else:
        tmpmol = molslist[0].copy()
    mpc = atomPerCell / len(tmpmol.sites)
    # check if mpc is an integer
    if mpc.is_integer():
        print('The number of molecules per cell is:', int(mpc))
    else:
        print('The number of molecules per cell is not an integer, exiting program...')
    return int(mpc)

def getEdgeFragmentsIndex(supercell, mollen, intermoldist, finegrid, bondDict, adjustment=1.0):
    print('Looking for the edge fragments index, might take up to an hour...')
    indexlist = list()
    # decide if cutoff need a different sign
    if max(supercell.frac_coords[:, 0]) < 0:
        sign = -1
    else:
        sign = 1
    # check in three dimensions
    for i in range(3):
        # cutoff = intermoldist * adjustment * sign / supercell.lattice.abc[i]
        cutoff = adjustment * sign * (1/finegrid[i])
        print('The cutoff is:', cutoff)
        # decreaselist tells the index for atoms within the cutoff range
        # search start from these atoms
        if sign < 0:
            decreaselist = np.where(np.logical_or(supercell.frac_coords[:, i] > cutoff, supercell.frac_coords[:, i] < (-1-cutoff)))[0]
        else:
            decreaselist = np.where(np.logical_or(supercell.frac_coords[:, i] < cutoff, supercell.frac_coords[:, i] > (1-cutoff)))[0]
        print('The length of decreaselist is:', len(decreaselist))
        # edgelist is a list of sites, it contains the sites belong to the edge fragments, not index
        # find out the index later
        edgeindex = list()
        # key function is: getSingleMol(supercell, middleSite, bondDict, middleSiteIndex)
        for cutindex in decreaselist:
            # if this index in the decreaselist is already in the cutoffindex
            # then don't waste time to find out the fragment for it
            if cutindex in edgeindex:
                pass
            else:
                # choose the first site index in the decreaselist as the starting middleSiteIndex
                # getSingleMol returns a dictionary, key is index, value is site
                fragment = getSingleMol(supercell, supercell.sites[cutindex], bondDict, cutindex)
                for key in fragment.keys():
                    edgeindex.append(key)
            # before remove the edge sites, make a copy
            print('The length of edge index is:', len(edgeindex))
        indexlist.append(edgeindex)
    return indexlist[0], indexlist[1], indexlist[2]

def getMoleculeLength(molslist):
    # make sure this list is not empty
    if len(molslist) == 0:
        raise Exception('The number of molecules in molslist is zero, \
            please check make sure your structure is correct')
    else:
        tmpmol = molslist[0].copy()
    # get the longest distance within one molecule
    dist = 0
    for i, sitex in enumerate(tmpmol.sites):
        for j, sitey in enumerate(tmpmol.sites):
            if i < j:
                tmpdist = np.linalg.norm(sitex.coords-sitey.coords)
                if tmpdist > dist:
                    dist = tmpdist
    # thus give the longest distance within one molecule
    return dist

def getBoxEdgeIndex(supercell, fineGrid, boxedgeDist):
    if max(supercell.frac_coords[:, 0]) < 0:
        sign = -1
    else:
        sign = 1
    # check three dimensions, boxedgeDist is the cutoff range
    indexlist = list()
    cutoff = boxedgeDist*sign
    print('The cutoff is:', cutoff)
    for i in range(3):
        if sign > 0:
            edgeIndex = np.where(np.logical_or(supercell.frac_coords[:, i] < cutoff, supercell.frac_coords[:, i] > (1-cutoff)))[0]
        elif sign < 0:
            edgeIndex = np.where(np.logical_or(supercell.frac_coords[:, i] > cutoff, supercell.frac_coords[:, i] < (-1-cutoff)))[0]
        indexlist.append(edgeIndex)
    return indexlist[0], indexlist[1], indexlist[2]

def getAllEdgeIndex(a, b, c):
    allindex = np.array([])
    allindex = np.append(allindex, a)
    allindex = np.append(allindex, b)
    allindex = np.append(allindex, c)
    uniqueIndex = np.unique(allindex)
    return uniqueIndex

def getIndexAroundHole(holeCoords, supercell, intermoldist, bondDict, moldistMul):
    # first find out the index of atoms within the sphere, with hole as the center
    # make sure include_index=True
    # so that the third returned value within a neighbor is the index
    sphereNeighbors = supercell.get_neighbors(supercell.sites[-1], intermoldist*moldistMul, include_index=True)
    # remove this He atom, which is actually the hole position
    supercell.remove_species(["He"])
    sphereFragmentIndex = list()
    for neighbor in sphereNeighbors:
        if neighbor[2] in sphereFragmentIndex:
            pass
        else:
            # tmpMol is a dict, with key:index, value:site
            tmpMol = getSingleMol(supercell, neighbor[0], bondDict, neighbor[2])
            for index in tmpMol.keys():
                sphereFragmentIndex.append(index)
    print('The length of complete fragments around the hole is:', len(sphereFragmentIndex))
    return sphereFragmentIndex
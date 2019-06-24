""" 
    utility.py provide all functions:

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
from itertools import filterfalse
import math
import json

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
        print()
        print('Reminder: Please make sure you type the correct fine grid')
        print()
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

def loadSingleMol(path):
    singleMol = Molecule.from_file(path+'/singlemolecule/singleMol.xyz')
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
    print()
    print('hole positions:')
    for key in holePositions.keys():
        print(key, ':', holePositions[key])
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
    print('charge site', chargeSite)
    print('hole position', holePosition)
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

def outputHolePositions(holeSites, path):
    supercellPath = path+'/supercell'
    # json does not allow numpy,int64
    # convert values into list of float
    filename = 'holePositions.json'
    if filename in os.listdir(supercellPath):
        print('There is one json file from previous calculation')
        print('The old file will be rewritten')
        decision = 'A+'
        while decision != 'Y' and decision != 'N':
            decision = input('Do you want to proceed? Y for \'yes\' and N for \'no\.')
            if decision == 'Y':
                pass
            elif decision == 'N':
                return 0
            else:
                print('Please type in either Y or N!!!')
    tmpdict = dict()
    for key in holeSites.keys():
        tmpdict[int(key)] = list(holeSites[key])
    with open(supercellPath+'/holePositions.json', 'w') as file:
        file.write(json.dumps(tmpdict, indent=4))

# this function will create a list of directories
# name them with the charge site index of the single molecule
# and place the plot_xct input files there
def createPlotxctInput(path, holeSites, plotxctinput):
    dbapath = os.path.join(path, 'dba')
    for key in holeSites.keys():
        plotxctpath = os.path.join(dbapath, str(key))
        os.system('mkdir '+plotxctpath)
        with open(os.path.join(plotxctpath, 'plotxct.inp'), 'w') as outfile:
            outfile.write(plotxctinput)
            outfile.write('hole_position   ')
            for value in holeSites[key]:
                outfile.write(str(value)+'   ')

def loadCubeCell(path):
    print('Loading supercell now, please wait...')
    namelist = os.listdir(path)
    asestruct = None
    for name in namelist:
        if name.endswith('cube'):
            asestruct = read(os.path.join(path, name))
            break
    if asestruct == None:
        print('Error!!!')
        print('There is no cube file under /data/supercell.')
        print('Please make sure the exciton wavefunction calculation output\
            is put under the folder \'supercell\'.')
        return 0
    else:
        pmgobj = AseAtomsAdaptor()
        pmgstruct = pmgobj.get_structure(asestruct)
        return pmgstruct

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
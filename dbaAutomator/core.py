""" 
    Author: Xingyu (Alfred) Liu 
    Email: xingyu.alfred.liu@gmail.com
    Description:
        core.py is the core of dbaAutomator, it provides two modules
        - automator: guide the dba process
        - checker: check convergence and CT character
"""

import os
from .functions import *
from .structio import *
from .ref import *
import time

# automator is the core class for dbaAutomator, functions are:
# 1. select the molecule in the middle of the supercell
# 2. get the hole positions for corresponding HOMO
# 3. calculate DBA according to previous procedures
class automator(object):

    def __init__(self, datapath, finegrid, filepath=None):
        self.path = datapath
        # check if there are necessary folders inside data folder
        checkDataFolder(self.path)
        if filepath != None:
            copyInput(filepath, self.path)
        self.fineGrid = finegrid
        self.bondCutoff = bondCutoff
        print('Now loading the unit cell information...')
        self.unitcell = loadUnitCell(self.path)
        tmpunitcell = self.unitcell.copy()
        self.supercell = getSuperCell(tmpunitcell, self.fineGrid)
        self.bondDict = getBondDict(self.unitcell, bondCutoff)
    
    def getmol(self, returnmol=False, outputmol=True):
        print('Now finding the central single molecule...')
        self.singleMol = getCentralSingleMol(self.supercell, self.bondDict)
        if outputmol:
            print()
            for key in self.singleMol.keys():
                print(key, self.singleMol[key])
            outputMolecule(self.singleMol, self.path)
        if returnmol:
            print('The chosen single molecule is:')
            for key in self.singleMol.keys():
                print(key, self.singleMol[key])
            return self.singleMol

    def getholes(self, returnholes=False, writeinput=True, chargeThreshold=0.01):
        self.chargeThreshold = chargeThreshold
        print('The charge threshold you chose for choosing the hole-placed site is:', "{:0.2f}".format(self.chargeThreshold), "%.")
        print('Now locate the hole positions...')
        print('Loading the output central single molecule...')
        self.centralmol = loadSingleMol(self.path)
        print('Loading cube and ACF.dat file...')
        self.chargeMatrix = loadChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
        print('Looking for hole positions...')
        tmpunitcell = self.unitcell.copy()
        self.holeSites = getHolePositions(self.chargeMatrix,self.centralmol,tmpunitcell,self.bondDict,self.chargeThreshold)
        if returnholes:
            print()
            print('The chosen holes are:')
            for key in self.holeSites.keys():
                print(key, ":", self.holeSites[key])
            return self.holeSites
        print('The hole positions are output under:', os.path.join(self.path, "supercell/holePositions.json"))
        outputHolePositions(self.holeSites, self.path)
        if writeinput:
            print('The input files for plotxct calculations are written under:', os.path.join(self.path, "dba"))
            createPlotxctInput(self.path, self.holeSites, self.fineGrid)

    def caldba(self, writeresult=True):
        print('Now calculate charge transfer character...')
        print('Loading the output central single molecule...')
        self.centralmol = loadSingleMol(self.path)
        print('Loading ACF.dat of single molecule HOMO...')
        self.smcharge = loadChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
        print('Loading bader results for each hole positions')
        dbapath = os.path.join(self.path, "dba")
        self.holeSites = loadHolePositions(self.path)
        # get the list of path for holes
        # need to check if they are all holes, delte those aren't
        holelist = os.listdir(dbapath)
        for hole in holelist[:]:
            if hole not in self.holeSites.keys():
                holelist.remove(hole)
        print('Checking if holes match with previous settings...')
        # check the other way
        for key in self.holeSites.keys():
            if key not in holelist:
                raise Exception('Hole position No.', key, 'is not in dba folder.')
        # the charge occupation for each hole will be stored at chargeshare dict
        chargeshare = dict()
        for hole in holelist:
            print()
            print('Loading ACF.dat at hole index:',hole)
            holepath = os.path.join(dbapath, hole)
            cubeSuperCell = loadCubeCell(holepath)
            chargematrix = loadChargeMatrix(cubeSuperCell, holepath)
            molIndex = getMoleculeIndex(self.centralmol, cubeSuperCell)
            chargeshare[hole] = getMolShare(chargematrix, molIndex)
            print('Charge occupation for molecule with hole', hole, 'is:', "{:0.2f}".format(chargeshare[hole]*100), "%.")
        print()
        print('The charge transfer character for each hole positions:')
        for hole in chargeshare.keys():
            print(hole, ":", "{:0.2f}".format((1-chargeshare[hole])*100), "%")
        print()
        print('Computing charge transfer character now...')
        self.holeindexlist = np.array(holelist).astype(int)
        tmpcharge = self.smcharge / np.sum(self.smcharge[self.holeindexlist])
        chargetransfer = 0
        for hole in self.holeindexlist:
            chargetransfer += chargeshare[str(hole)] * tmpcharge[hole][4]
        print('The total charge transfer character is:', "{:0.2f}".format((1-chargetransfer)*100), "%.")
        # at last, decide if want to write the result into a file
        if writeresult:
            writedbaResult(self.path, chargeshare, chargetransfer)

# checker is a supporting class, it can check the convergence for plotxct calculation
# it can also calculate CT character for individual folders
class checker(object):

    # the path has to be the absolute path to where dba runs
    def __init__(self, path, finegrid):
        self.path = path
        self.fineGrid = finegrid
        checkList = []
        self.checklist = getXctPath(self.path, checkList)
        if len(self.checklist) == 0:
            raise Exception('Warning!!! No suitable files found in folder:', self.path, "\n")
    
    def prep(self):
        # choose a supercell cell and get the inter molecular distance
        # first get all complete single molecules
        tmpsupercell = loadCubeCell(self.checklist[0])
        supercell = tmpsupercell.copy()
        self.bondDict = getBondDict(tmpsupercell, bondCutoff)
        print('Looking for the primitive cell...')
        self.primitive = getPrimitiveCell(tmpsupercell)
        print('Looking for all fragments within constructed supercell...')
        tmpstruct = self.primitive.copy()
        tmpstruct = getSuperCell(tmpstruct, [2, 2, 2])
        self.molslist = getAllMols(tmpstruct, self.bondDict)
        self.intermoldist = getInterMolLen(self.molslist)
        self.mollen = getMoleculeLength(self.molslist)
        print('The longest distance in one molecule is:', "{:0.2f}".format(self.mollen))
        self.mpc = getMPC(supercell, self.fineGrid, self.molslist)
        print('The closest distance between molecular center of masses is:', "{:0.2f}".format(self.intermoldist))

    def checkconv(self, convThreshold=0.05, convDist=1.0):
        start_time = time.time()
        print('Checking the convergence for founded exciton wavefunction calculations...')
        # need to get the index for the cube file edge fragments
        supercell = loadCubeCell(os.path.join(self.checklist[0]))
        # print('Getting the index for edge fragments')
        # self.edgeAindex, self.edgeBindex, self.edgeCindex = getEdgeFragmentsIndex(supercell, self.mollen,\
        #     self.intermoldist, self.fineGrid, self.bondDict, adjustment=convDist)
        print('time cost is:', time.time()-start_time)
        for name in self.checklist:
            os.chdir(name)
            print()
            print('Now check folder:', name)
            # print('Loading ACF.dat...')
            print('Loading ACF.dat and cube file...')
            self.boxedgeAindex, self.boxedgeBindex, self.boxedgeCindex = getBoxEdgeIndex(supercell, self.fineGrid, boxedgeDist=convDist)
            chargematrix = loadChargeMatrix(supercell, name)
            # chargedira = getChargeShare(self.edgeAindex, chargematrix)
            # chargedirb = getChargeShare(self.edgeBindex, chargematrix)
            # chargedirc = getChargeShare(self.edgeCindex, chargematrix)
            chargedira = getChargeShare(self.boxedgeAindex, chargematrix)
            chargedirb = getChargeShare(self.boxedgeBindex, chargematrix)
            chargedirc = getChargeShare(self.boxedgeCindex, chargematrix)
            printChargeShare(chargedira, chargedirb, chargedirc, convThreshold)
            self.boxEdgeAll = getAllEdgeIndex(self.boxedgeAindex, self.boxedgeBindex, self.boxedgeCindex)
            print(self.boxEdgeAll)
            chargeAllEdge = getChargeShare(self.boxEdgeAll.astype(int), chargematrix)
            print('The total charge share for all edge sites is:', "{:0.2f}".format(chargeAllEdge*100), "%")
            os.chdir('../')

    def calct(self, filepath):
        self.unitcell = loadUnitCell(filepath)
        self.bondDict = getBondDict(self.unitcell, bondCutoff)
        for name in self.checklist:
            os.chdir(name)
            print()
            print('Now calculating charge transfer character in folder:', name)
            supercell = loadCubeCell(name)
            chargematrix = loadChargeMatrix(supercell, name)
            holePosition = loadPlotxct(name)
            print('The hole position in the input file is:', holePosition)
            tmpunitcell = self.unitcell.copy()
            tmpunitcell.append('He', holePosition)
            holeCartesianCoords = tmpunitcell.sites[-1].coords
            print('The Cartesian coordinates for this hole position is:', holeCartesianCoords)
            # get teh fractional coords for hole
            tmpcube = supercell.copy()
            tmpcube.append('He', holeCartesianCoords, coords_are_cartesian=True)
            holeFracCoords = tmpcube.sites[-1].frac_coords
            singleMol = getCentralSingleMol(supercell, self.bondDict, middle=holeFracCoords)
            # need to construct a Molecule object to be passed into getMoleculeIndex
            centralmol = Molecule([], [])
            for siteindex in singleMol.keys():
                centralmol.append(str(singleMol[siteindex].specie), singleMol[siteindex].coords)
            molIndex = getMoleculeIndex(centralmol, supercell)
            # calculate the Frenkel character for this hole position
            frenkel = getMolShare(chargematrix, molIndex)
            print('The charge transfer character for this hole position is:', "{:0.2f}".format((1-frenkel)*100), "%.")
            os.chdir("../")
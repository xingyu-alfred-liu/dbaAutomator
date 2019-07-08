import os
from functions import *
from structio import *
from ref import *

class automator(object):
    def __init__(self, path, finegrid, chargeThreshold=0.01):
        self.path = path
        self.fineGrid = finegrid
        self.bondCutoff = bondCutoff
        self.chargeThreshold = chargeThreshold
        print('Now loading the unit cell information...')
        self.unitcell = loadUnitCell(self.path)
        self.supercell = getSuperCell(self.unitcell, self.fineGrid)
        self.bondDict = getBondDict(self.unitcell, bondCutoff)
    
    def getcentralmol(self, returnmol=False, outputmol=True):
        print('Now finding the central single molecule...')
        self.singleMol = getCentralSingleMol(self.supercell, self.bondDict)
        if outputmol:
            outputMolecule(self.singleMol, self.path)
        if returnmol:
            print('The chosen single molecule is:')
            print(self.singleMol)
            return self.singleMol

    def getholes(self, returnholes=False, outputholes=True, writeinput=True):
        print('Now locate the hole positions...')
        print('Loading the output central single molecule...')
        self.centralmol = loadSingleMol(self.path)
        print('Loading ACF.dat file...')
        self.chargeMatrix = loadChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
        print('Looking for hole positions...')
        self.holeSites = getHolePositions(self.chargeMatrix,self.centralmol, \
                                          self.unitcell, self.bondDict, self.chargeThreshold)
        if returnholes:
            print()
            print('The chosen holes are:')
            for key in self.holeSites.keys():
                print(key, ":", self.holeSites[key])
            return self.holeSites
        if outputholes:
            print('The hole positions are output under:', os.path.join(self.path, "supercell/holePositions.json"))
            outputHolePositions(self.holeSites, self.path)
        if writeinput:
            print('The input files for plotxct calculations are written under:', os.path.join(self.path, "dba"))
            createPlotxctInput(self.path, self.holeSites, self.fineGrid)

    def caldba(self):
        print('Now calculate charge transfer character...')
        print('Loading the output central single molecule...')
        self.centralmol = loadSingleMol(self.path)
        print('Loading ACF.dat of single molecule HOMO...')
        self.smcharge = loadChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
        print('Loading bader results for each hole positions')
        dbapath = os.path.join(self.path, "dba")
        self.holeSites = loadHolePositions(self.path)
        holelist = os.listdir(dbapath)
        print('Checking if holes match with previous settings...')
        for key in self.holeSites.keys():
            if key not in holelist:
                raise Exception('Hole position No.', key, 'is not in dba folder.')
        # the charge occupation for each hole will be stored at chargeshare dict
        chargeshare = dict()
        for hole in holelist:
            print('Loading supercell and ACF.dat at hole index:',hole)
            holepath = os.path.join(dbapath, hole)
            cubeSuperCell = loadCubeCell(holepath)
            chargematrix = loadChargeMatrix(cubeSuperCell, holepath)
            molIndex = getMoleculeIndex(self.centralmol, cubeSuperCell)
            chargeshare[hole] = getMolShare(chargematrix, molIndex)
            print('Charge occupation for hole index', hole, 'is:', chargeshare[hole])
        print('The charge transfer character for each hole positions:')
        for hole in chargeshare.keys():
            print(hole, ":", chargeshare[hole])
        print('Computing charge transfer character now...')
        holeindexlist = np.array(holelist).astype(int)
        self.smcharge = self.smcharge / (np.sum(self.smcharge(holeindexlist), axis=0)[4])
        chargetransfer = 0
        for hole in holeindexlist:
            chargetransfer += chargeshare[str(hole)] * self.smcharge[hole]
        print('The total charge transfer character is:', "{:0.2f}".format(chargetransfer*100), "%.")

class checker(object):
    # the path has to be the absolute path to where dba runs
    def __init__(self, path):
        self.path = path
        checkList = []
        self.checklist = getXctPath(self.path, checkList)
        if len(self.checklist) == 0:
            raise Exception('Warning!!! No suitable files found in folder:', self.path, "\n")
        # choose a supercell cell and get the inter molecular distance
        # first get all complete single molecules
        self.tmpsupercell = loadCubeCell(self.checklist[0])
        self.bondDict = getBondDict(self.tmpsupercell, bondCutoff)
        print('Looking for the primitive cell...')
        self.primitive = getPrimitiveCell(self.tmpsupercell)
        print('Looking for all fragments within constructed supercell...')
        tmpstruct = self.primitive.copy()
        tmpstruct = getSuperCell(tmpstruct, [2, 2, 2])
        self.molslist = getAllMols(tmpstruct, self.bondDict)
        self.convrange = getInterMolLen(self.molslist)
        print('The closest distance between center of masses is:', "{:0.2f}".format(self.convrange))
    def checkconv(self, convThreshold=0.1):
        for name in self.checklist:
            os.chdir(name)
            print('Now check folder:', name)
            chargematrix = loadChargeMatrix(self.tmpsupercell, name)
            moldira, moldirb, moldirc = getAtomIndex(self.tmpsupercell, self.convrange)
            chargedira = getChargeShare(moldira, chargematrix)
            chargedirb = getChargeShare(moldirb, chargematrix)
            chargedirc = getChargeShare(moldirc, chargematrix)
            printChargeShare(chargedira, chargedirb, chargedirc, convThreshold)
            os.chdir('../')
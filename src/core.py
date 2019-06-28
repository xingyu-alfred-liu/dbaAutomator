import os
from functions import *
from ref import *

class automator(object):
    def __init__(self, path, fineGrid, chargeThreshold=0.01):
        self.path = path
        self.fineGrid = fineGrid
        self.bondCutoff = bondCutoff
        self.plotxctinput = plotxctinput
        self.chargeThreshold = chargeThreshold
        print('Now loading the unit cell information...')
        self.unitcell = loadUnitCell(self.path)
        self.supercell = getSuperCell(self.path, self.unitcell, self.fineGrid)
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
        self.chargeMatrix = getChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
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
            createPlotxctInput(self.path, self.holeSites, self.plotxctinput)

    # the path has to be the path to a folder where one cube file and corresponding 
    def checkconvergence(self, path):
        print("Checking convergence now...")
        filelist = os.listdir(path)
        for file in filelist:
            if file.endswith('.cube'):
                cubefilename = file
                break
        cubeSupercell_check = loadCubeCell(os.path.join(path, cubefilename))
        chargeMatrix_check = getChargeMatrix(cubeSupercell_check, os.path.join(path, path))

    def caldba(self):
        print('Now calculate charge transfer character...')
        print('Loading the output central single molecule...')
        self.centralmol = loadSingleMol(self.path)
        print('Loading ACF.dat of single molecule HOMO...')
        self.smcharge = getChargeMatrix(self.centralmol, os.path.join(self.path, 'singlemolecule'))
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
            chargematrix = getChargeMatrix(cubeSuperCell, holepath)
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
        print('The total charge transfer character is:', "{:0.2f}".format(chargetransfer*100), "%")
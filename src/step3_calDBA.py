from utility import *
from files import *
import os

class dbaautomator(object):
    def __init__(self, path):
        self.path = path
    def loadSingleMol(self):
        self.singleMol = loadSingleMol(self.path)
    def loadCubeCell(self):
        self.cubeSuperCell = loadCubeCell(os.path.join(self.path, 'supercell'))
    def getMoleculeIndex(self):
        self.molIndex = getMoleculeIndex(self.singleMol, self.cubeSuperCell)
    def getChargeMatrix(self):
        self.chargeMatrix = getChargeMatrix(self.cubeSuperCell, os.path.join(self.path, 'supercell'))
    def getMolShare(self):
        self.molShare = getMolShare(self.chargeMatrix, self.molIndex)

if __name__ == "__main__":
    singleMol = loadSingleMol(dataDir)
    cubeSuperCell = loadCubeCell(os.path.join(dataDir, 'supercell'))
    molIndex = getMoleculeIndex(singleMol, cubeSuperCell)
    chargeMatrix = getChargeMatrix(cubeSuperCell, os.path.join(dataDir, 'supercell'))
    molShare = getMolShare(chargeMatrix, molIndex)
    print('The charge share for this site is:', "{:0.2f}".format(molShare*100), "%")
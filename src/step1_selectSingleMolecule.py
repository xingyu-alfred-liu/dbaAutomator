"""
    The first step is to locate the single molecule cloest to the middle of the 
    constructed supercell.
    The reason for this is to help locate the hole positioned molecule, and help
    calculate the Frenkel character. 

"""

from utility import *
from files import *
from bondLengthRef import bondCutoff

class singleMolSelector(object):
    def __init__(self, path):
        self.path = path
    def loadUnitCell(self):
        print('Step 1: Find the single molecule cloest to the middle in super cell')
        self.unitcell = loadUnitCell(self.path)
        print('Done with loading unit cell')
    def getSuperCell(self, fineGrid):
        self.supercell = getSuperCell(self.path, self.unitcell, fineGrid)
    def getbondDict(self, bondCutoff):
        self.bondDict = getBondDict(self.supercell, bondCutoff)
    def getCentralSingleMol(self):
        self.singleMol = getCentralSingleMol(self.supercell, self.bondDict)
    def outputMolecule(self):
        outputMolecule(self.singleMol, self.path)

def step1(dataDir, fineGrid, bondCutoff):
    # print('Step 1: Find the single molecule cloest to the middle in super cell')
    # unitcell = loadUnitCell(dataDir)
    # supercell = getSuperCell(dataDir, unitcell, fineGrid)
    # bondDict = getBondDict(supercell, bondCutoff)
    # singleMol = getCentralSingleMol(supercell, bondDict)
    # outputMolecule(singleMol, dataDir)
    smSelector = singleMolSelector(dataDir)
    smSelector.loadUnitCell()
    smSelector.getSuperCell(fineGrid)
    smSelector.getbondDict(bondCutoff)
    smSelector.getCentralSingleMol()
    smSelector.outputMolecule()

if __name__ == "__main__":
    step1(dataDir, fineGrid, bondCutoff)
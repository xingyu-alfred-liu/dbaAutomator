from utility import *
from bondLengthRef import bondCutoff

class holeLocator(object):
    def __init__(self, path):
        self.path = path
    def loadSingleMol(self):
        self.singleMol = loadSingleMol(self.path)
    def getChargeMatrix(self):
        self.chargeMatrix = getChargeMatrix(self.singleMol, os.path.join(self.path, 'singlemolecule'))
    def loadUnitCell(self):
        self.unitcell = loadUnitCell(self.path)
    def getBondDict(self, bondCutoff):
        self.bondDict = getBondDict(self.unitcell, bondCutoff)
    def getHolePositions(self, chargeThreshold):
        self.holeSites = getHolePositions(self.chargeMatrix,\
            self.singleMol, self.unitcell, self.bondDict, chargeThreshold)
    def outputHolePositions(self):
        outputHolePositions(self.holeSites, self.path)
    def createPlotxctInput(self, plotxctinput):
        createPlotxctInput(self.path, self.holeSites, plotxctinput)

if __name__ == "__main__":
    # singleMol = loadSingleMol(dataDir)
    # chargeMatrix = getChargeMatrix(singleMol, dataDir)
    # unitcell = loadUnitCell(dataDir)
    # bondDict = getBondDict(unitcell, bondCutoff)
    # holeSites = getHolePositions(chargeMatrix, singleMol, unitcell, bondDict, chargeThreshold)
    holelocator = holeLocator(dataDir)
    holelocator.loadSingleMol()
    holelocator.getChargeMatrix()
    holelocator.loadUnitCell()
    holelocator.getBondDict(bondCutoff)
    holelocator.getHolePositions(chargeThreshold)
    holelocator.outputHolePositions()
    holelocator.createPlotxctInput(plotxctinput)
from files import *
from utility import *
from bondLengthRef import bondCutoff

if __name__ == "__main__":
    singleMol = loadSingleMol(dataDir)
    chargeMatrix = getChargeMatrix(singleMol, dataDir)
    unitcell = loadUnitCell(dataDir)
    bondDict = getBondDict(unitcell, bondCutoff)
    holeSites = getHolePositions(chargeMatrix, singleMol, unitcell, bondDict, chargeThreshold)

"""
    - ckecking if the excition wfn calculation is converged is
    checking if 1/8 of the supercell occupy more than 90% of the 
    excited electron
"""
from filepath import *
from calMolShare import getSuperCell
from locateHole import getChargeMatrix
from numpy import np

def getchargeoccupation(supercellPath):
    cubeSuperCell = getSuperCell(supercellPath)
    chargeMatrix = getChargeMatrix(cubeSuperCell, supercellPath)
    atomRange = [0.25, 0.75]
    return None

if __name__ == "__main__":
    chargeOccupation = getchargeoccupation(supercellPath)

"""
    - ckecking if the excition wfn calculation is converged is
    checking if 1/8 of the supercell occupy more than 90% of the 
    excited electron
"""
from filepath import *
from calMolShare import getSuperCell
from locateHole import getChargeMatrix
import numpy as np

def getChargeOccupation(cubeSuperCell, chargeMatrix):
    atomRange = [0.25, 0.75]
    atomSiteIndex = list()
    for i, site in enumerate(cubeSuperCell.sites):
        if len(np.where(site.frac_coords>atomRange[0])[0])==3 and len(np.where(site.frac_coords<atomRange[1])[0])==3:
            atomSiteIndex += [i]
    chargeSum = np.sum(chargeMatrix[atomSiteIndex], axis=0)[4]
    return chargeSum

def outputConvergenceCondition(chargeOccupation):
    print("1/8 of the total supercell occupies:", "{:0.2f}".format(chargeOccupation*100), "%", "charge.")
    if chargeOccupation < 0.9:
        print("Warning!!!")
        print("This exciton wavefunction might not be converged. Please check your calculation.")
    else:
        print("Your calculation is converged. Please proceed.")

if __name__ == "__main__":
    cubeSuperCell = getSuperCell(supercellPath)
    chargeMatrix = getChargeMatrix(cubeSuperCell, supercellPath)
    chargeOccupation = getChargeOccupation(cubeSuperCell, chargeMatrix)
    outputConvergenceCondition(chargeOccupation)
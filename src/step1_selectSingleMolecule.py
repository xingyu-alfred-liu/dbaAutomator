"""
    The first step is to locate the single molecule cloest to the middle of the 
    constructed supercell.
    The reason for this is to help locate the hole positioned molecule, and help
    calculate the Frenkel character. 

"""

from utility import *
from files import *
from bondLengthRef import bondCutoff

if __name__ == "__main__":
    print('Step 1: Find the single molecule cloest to the middle in super cell')
    unitcell = loadUnitCell(dataDir)
    print(unitcell)
    supercell = getSuperCell(dataDir, unitcell, fineGrid)
    bondDict = getBondDict(supercell, bondCutoff)
    singleMol = getCentralSingleMol(supercell, bondDict)
    outputMolecule(singleMol, dataDir)
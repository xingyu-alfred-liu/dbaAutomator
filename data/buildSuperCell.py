from ase.io import read,write
from ase.build import make_supercell

inputfile = 'in'
structure = read(inputfile)
supercell = make_supercell(structure, [[8, 0, 0], [0, 8, 0], [0, 0, 2]])
print(supercell)

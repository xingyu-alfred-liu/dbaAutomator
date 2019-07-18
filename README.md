# dbaAutomator
## Introduction

dbaAutomator automates double-bader analysis (DBA). It interfaces with exciton wavefunction (ext wfn) calculations output `cube` files performed by [BerkeleyGW](https://berkeleygw.org/), orbital files computed with most DFT codes such as [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/) and [Quantum ESPRESSO](https://www.quantum-espresso.org/), [Bader Charge Analysis](http://theory.cm.utexas.edu/henkelman/code/bader/) output `ACF.dat`.

**_Note: we are still working on extending the functionalities of dbaAutomator, such as provide bader analysis binaries, compute HOMO for selected single molecule, etc._**

#### Requirements
  * numpy
  * ase
  * pymatgen  

#### Modules
  * **automator** - provide dba workflow guidance, including output the hole-placed single molecule, output hole positions, generate inputs for `plotxct.x`, and compute separate and total charge transfer character.
  * **checker** - compute intermolecular distance, check convergence of excitation calculation, compute charge transfer character with give bader analysis output.

## Installation

### Install from source

#### Clone the repository
`git clone https://github.com/BLABABA/dbaAutomator.git`  
`cd dbaAutomator`

#### Install requirements
`pip install -r requirements.txt`  
  
Notice that `pymatgen` might not be installed with `pip`, the other choice is using `conda`  

`conda install --channel matsci pymatgen`  
  
In case you are not familiar with `pymatgen`, here is the [link](https://pymatgen.org/)

`numpy` and `ase` can be installed with `pip`

#### Install dbaAutomator
`python setup.py install`  
`cd ..`  
  
**You are ready to take off ;)**

## Tutorial

### Guide the procedure of DBA

#### Define an automator object
Import `dbaAutomator.core.automator` and `os`  

    >>> import os
    >>> from dbaAutomator.core import automator

*__class__* `automator(datapath, finegrid, filepath=None, chargeThreshold=0.01)`  
Description: guide the dba process.

**Parameters**  

`datapath`: **_string_**  
The absolute path where you want to put your data in. We suggest every users create a new directory, such as:  

    >>> os.mkdir('data')
    >>> os.chdir('data')
    >>> os.path.abspath('.')
    '/TheBestUser/PATH/TO/data'

`finegrid`: **_list_**  
The fine grid defined in mean-field calculation to obtain the fine grid wavefunction `WFN`. It can be found in file `kgrid.in`. Here we define `[a, b, c]` as a list of arbitrary positive integers.

    >>> finegrid = [a, b, c]
`filepath`: **_string_**  
Defulat: `None`  
The absolute path to the input file which defines the unit cell. Here we use the `in` of Quantum ESPRESSO as an example. If not set, input will be only searched under `datapath`.

    >>> filepath = '/TheBestUser/PATH/TO/in'
  
Now one can define the automator object. The initial setup will lead to creation of necessary directories, duplication of input file to directory `unitcell`, and construction of supercell. 

    >>>> dba = automator(path=datadir, finegrid=finegrid, file=filepath)
    Directory dba created.
    Directory singlemolecule created.
    Directory supercell created.
    Directory unitcell created.
    Copy unitcell file to /Users/alfredliu/tmp/unitcell .
    Now loading the unit cell information...
    Please make sure this is the file you want to load: /Users/alfredliu/tmp/unitcell/in
    Reminder: Please make sure you type the correct fine grid

#### Step 1: Find the single molecule

*__function__* `getmol(returnmol=False, outputmol=True)`  
Description: select the single molecule in the middle of the constructed supercell.

**Parameters**  
`returnmol`: **_bool_**  
Default: `False`  
Return a `pymatgen.Molecule` object of this selected molecule.  

`outputmol`: **_bool_**  
Default: `True`  
The single molecule file will be written in the singlemol folder. 

    >>> dba.getmol(returnmol=False, outputmol=True)
    Now finding the central single molecule...
    The site closest to the middle is [  -5.38261049   -8.45632748 -152.19234847] H
    The corresponding site index is 6831
    Number of atoms found in this molecule: 1
    Number of atoms found in this molecule: 2
    Number of atoms found in this molecule: 4
    ......
    Number of atoms found in this molecule: 122
    Number of atoms found in this molecule: 133
    Number of atoms found in this molecule: 136
    
    6831 [  -5.38261049   -8.45632748 -152.19234847] H
    34479 [  -4.33330991   -8.20105926 -152.05995518] C
    28335 [  -3.95927738   -6.89532748 -151.71176732] C
    ......
    7623 [  -6.72688503    1.33136346 -160.59428942] H
    14703 [  -5.16135803   -5.60808548 -139.17673664] H
    16239 [  -5.71045515   -3.38168353 -139.97871313] H
    The single molecule structure is saved under '/TheBestUser/PATH/TO/data/singlemolecule/singleMol.xyz'

#### Step 2: Locate the hole positions

*__function__* `getholes(returnholes=False, writeinput=True, chargeThreshold=0.01)`  
Description: locate the hole positions and write the input files for ext wfn calculations and write the hole positions into a json file. 

**Parameters**  

`returnholes`: **_bool_**  
Default: `True`  
Decide if return a dictionary, with keys as the hole-placed atom index, and values are the fractional coords of holes with respect to the unit cell.  
**Notice**: this file will be loaded in dba calculation, so please keep it as where it is.  

`writeinput`: **_bool_**  
Default: `True`  
The input files for ext wfn calculations will be generated under:  
`/TheBestUser/PATH/TO/data/dba`  
Each folder within this path represents a hole position.  
**Notice**: Please keep it as default `True`, because these folders will be loaded in dba calculation. 

`chargeThreshold`: **_float_**  
Default: `0.01`  
Define the charge threshold of possible hole-placed atom sites. In other words, if one atom occupies charge more than this threshold, it will be considered as possible position to put hole near to. 

    >>> dba.getholes(returnholes=False, writeinput=True)
    Now locate the hole positions...
    Loading the output central single molecule...
    Loading cube and ACF.dat file...
    Looking for hole positions...
    The hole positions are output under: /TheBestUser/PATH/TO/data/supercell/holePositions.json
    The input files for plotxct calculations are written under: /TheBestUser/PATH/TO/data/dba

#### Step 3: Compute charge transfer character
*__function__* `caldba(writeresult=True)`  
Description: compute the charge transfer character for each hole position, and the total charge transfer character.  

**Parameters**  

`writeresult`: **_bool_**  
Default: `True`  
The charge transfer character information will be output under: `/TheBestUser/PATH/TO/data/supercell`

    >>> dba.caldba(writeresult=False)
    Now calculate charge transfer character...
    Loading the output central single molecule...
    Loading ACF.dat of single molecule HOMO...
    Loading bader results for each hole positions
    Checking if holes match with previous settings...
    
    Loading ACF.dat at hole index: 66
    Loading supercell, this process might take minutes, please wait...
    Charge occupation for hole index 66 is: 0.60 %.
    ......
    Loading ACF.dat at hole index: 87
    Loading supercell, this process might take minutes, please wait...
    Charge occupation for hole index 87 is: 0.49 %.

    The charge transfer character for each hole positions:
    66 : 99.40 %
    ......
    87 : 99.51 %
    
    Computing charge transfer character now...
    The total charge transfer character is: 99.30 %.


### Check convergence, charge transfer

#### Define the checker object
Import `dbaAutomator.core.checker`  

*__class__* `checker(path)`  
Description: used to check convergence of ext wfn calculation and calculate charge transfer character for an individual result. 

**Parameters**  

`path`: **_string_**  
Feed in the path where you want to run checker. Check will run recurrently and find out all folder with `cube` file and corresponding `ACF.dat`, appended as part of a checklist. 

    >>> from dbaAutomator.core import checker
    >>> checkpath = '/TheBestUser/PATH/TO/data/dba'
    >>> dbachecker = checker(checkpath)

#### Compute the intermolecular distance

*__function__* `caldist()`  
Description: load in the cube supercell and find all complete single molecules, and find out the smallest intermolecular distance. 

    >>> dbachecker.caldist()
    Loading supercell, this process might take minutes, please wait...
    Looking for the primitive cell...
    Looking for all fragments within constructed supercell...
    Reminder: Please make sure you type in the correct fine grid.
    length of supercell: 2176
    The site closest to the middle is [  3.99770256  19.97457752 -15.96009805] C
    The corresponding site index is 450
    Number of atoms found in this molecule: 1
    Number of atoms found in this molecule: 4
    ...
    Number of atoms found in this molecule: 132
    Number of atoms found in this molecule: 136
    length of supercell: 2040
    The site closest to the middle is [  3.12439252  21.83230651 -18.94876799] H
    The corresponding site index is 122
    Number of atoms found in this molecule: 1
    Number of atoms found in this molecule: 2
    Number of atoms found in this molecule: 4
    ...
    length of supercell: 136
    The site closest to the middle is [  1.21158426   6.77991156 -26.88168128] H
    The corresponding site index is 114
    Number of atoms found in this molecule: 1
    Number of atoms found in this molecule: 2
    ...
    Number of atoms found in this molecule: 134
    Number of atoms found in this molecule: 136
    The closest distance between center of masses is: 8.21

#### Check convergence of ext wfn calculation

*__function__* `checkconv(convThreshold=0.05)`  
Description: load in the cube supercell and find all complete single molecules, and find out the smallest intermolecular distance. 

**Parameters**

`convThreshold `: **_float_**  
Default: `0.05`
Define the convergence threshold of charge occupied by edge atoms. 

#### Compute charge transfer for specified conditions

*__function__* `calct(filepath)`  
Description:  compute the charge transfer character for ext wfn calculation output defined in `check` class. 

**Parameters**

`filepath`: **_string_**  
Defines the path to the mean field calculation input, which provides the unit cell structure. 

    >>> dbachecker.calct('/TheBestUser/PATH/TO/in')
    
    Now calculating charge transfer character in folder: /CHECK/LIST/66
    Loading supercell, this process might take minutes, please wait...
    The hole position in the input file is: [3.141523121238201, 2.0205509312877465, 2.934486265836001]
    The Cartesian coordinates for this hole position is: [15.26745904 42.93751551 46.18576196]
    The site closest to the middle is [15.75578418 43.25327987 45.63637017] C
    The corresponding site index is 23626
    Number of atoms found in this molecule: 1
    ...
    Number of atoms found in this molecule: 134
    Number of atoms found in this molecule: 136
    The charge transfer character for this hole position is: 99.40 %.
    
    Now calculating charge transfer character in folder: /CHECK/LIST/15
    Loading supercell, this process might take minutes, please wait...
    ...
    

## External Resources:  
  * Cube File explanation: [http://paulbourke.net/dataformats/cube/](http://paulbourke.net/dataformats/cube/)
  * Bader Charge Analysis: [http://theory.cm.utexas.edu/henkelman/code/bader](http://theory.cm.utexas.edu/henkelman/code/bader)

## Notes:
DO NOT rerun single molecule selection after you get the HOMO, because the order of this single molecule might change.  
Do rerun single molecule selection if you decide to change the supercell size.

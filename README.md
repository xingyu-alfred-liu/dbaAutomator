# dbaAutomator
## Introduction

dbaAutomator automates double-bader analysis (DBA). It interfaces with exciton wavefunction (ext wfn) calculations output `cube` files performed by [BerkeleyGW]
(https://berkeleygw.org/), orbital files computed with most DFT codes such as [FHI-aims](https://aimsclub.fhi-berlin.mpg.de/) and [Quantum ESPRESSO](https://www.quantum-espresso.org/), [Bader Charge Analysis](http://theory.cm.utexas.edu/henkelman/code/bader/) output `ACF.dat`.

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
  
**You are ready to take off.**

## Tutorial

### Guide the procedure of DBA

#### Define an automator object
Import `dbaAutomator.core.automator` and `os`  

    >>> import os
    >>> from dbaAutomator.core import automator

Define `datapath, finegrid, filepath`  

`datapath`: **_string_**  
The absolute path where you want to put your data in. We suggest every users create a new directory, such as:  

    >>> os.mkdir data
    >>> os.chdir('data')
    >>> os.path.abspath('.')
    '/TheBestUser/PATH/TO/data'

`finegrid`: **_list_**  
The fine grid defined in mean-field calculation to obtain the fine grid wavefunction `WFN`. It can be found in file `kgrid.in`. Here we define `[a, b, c]` as a list of arbitrary positive integers.

    >>> finegrid = [a, b, c]
`filepath`: **_string_**  
The absolute path to the input file which defines the unit cell. Here we use the `in` of Quantum ESPRESSO as an example.

    >>> filepath = '/TheBestUser/PATH/TO/in'
Define the automator object. The initial setup will lead to creation of necessary directories, duplication of input file to directory `unitcell`, and construction of supercell. 

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
Look for the single molecule located in the middle of the constructed supercell.  

**Arguments**  

`returnmol=True`  Default: `False`  
`dba.getmol()` will return a `pymatgen.Molecule` object of this selected molecule.  

`outputmol=True`  Default: `True`  
The single molecule file will be written in the singlemol folder. 

    >>> dba.getmol(returnmol=False, outputmol=True)
    Now finding the central single molecule...
    The site closest to the middle is [  -5.38261049   -8.45632748 -152.19234847] H
    The corresponding site index is 6831
    Number of atoms found in this molecule: 1
    Number of atoms found in this molecule: 2
    Number of atoms found in this molecule: 4
    Number of atoms found in this molecule: 8
    Number of atoms found in this molecule: 13
    Number of atoms found in this molecule: 18
    Number of atoms found in this molecule: 25
    Number of atoms found in this molecule: 32
    Number of atoms found in this molecule: 40
    Number of atoms found in this molecule: 54
    Number of atoms found in this molecule: 75
    Number of atoms found in this molecule: 101
    Number of atoms found in this molecule: 122
    Number of atoms found in this molecule: 133
    Number of atoms found in this molecule: 136
    
    6831 [  -5.38261049   -8.45632748 -152.19234847] H
    34479 [  -4.33330991   -8.20105926 -152.05995518] C
    28335 [  -3.95927738   -6.89532748 -151.71176732] C
    51573 [  -3.37141158   -9.19015407 -152.26717613] C
    45615 [  -5.13834955   -5.99033753 -140.19603869] C
    
    ......
    
    47343 [  -5.73454952   -3.745924   -141.00478986] C
    17007 [  -6.44231229   -1.89253381 -141.82349253] H
    4336 [ -10.12117455   -0.6406022  -145.63013462] H
    19503 [  -4.48778153   -9.85853304 -141.15639721] H
    18927 [  -4.80246479   -8.01226689 -139.60897644] H
    7623 [  -6.72688503    1.33136346 -160.59428942] H
    14703 [  -5.16135803   -5.60808548 -139.17673664] H
    16239 [  -5.71045515   -3.38168353 -139.97871313] H
    The single molecule structure is saved under '/TheBestUser/PATH/TO/data/singlemolecule/singleMol.xyz'

#### Step 2: Locate the hole positions
`dba.getholes()` can locate the hole positions and write the input files for ext wfn calculations under the path:  
`/TheBestUser/PATH/TO/data/dba`.  

**Arguments**  

`returnholes=True`  Default: `True`  

`dba.getholes()` will return a dictionary, with keys as the hole-placed atom index, and values are the fractional coords of holes with respect to the unit cell. A json file will be written under:  
`/TheBestUser/PATH/TO/data/supercell/holePositions.json`  
which contains the hole position information.  
Notice this file will be loaded in dba calculation, so please keep it as where it is.  

`writeinput=True`   Default: `True`  
The input files for ext wfn calculations will be generated under:  
`/TheBestUser/PATH/TO/data/dba`  
Each folder within this path represents a hole position. Please keep it as default `True`, because these folders will be loaded in dba calculation. 

    >>> dba.getholes(returnholes=False, outputholes=True, writeinput=True)
    Now locate the hole positions...
    Loading the output central single molecule...
    Loading cube and ACF.dat file...
    Looking for hole positions...
    The hole positions are output under: /TheBestUser/PATH/TO/data/supercell/holePositions.json
    The input files for plotxct calculations are written under: /TheBestUser/PATH/TO/data/dba

#### Step 3: Compute charge transfer character
The last step is using `dba.caldba()` to compute the charge transfer character for each hole position, and the total charge transfer character.  

**Arguments**  
`writeresult=True`  Default: `True`  
The charge transfer character inforlation will be output under: `/TheBestUser/PATH/TO/data/supercell`

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

#### Compute the intermolecular distance

#### Check convergence of ext wfn calculation

#### Compute charge transfer for specified conditions

## External Resources:  
  * Cube File explanation: [http://paulbourke.net/dataformats/cube/](http://paulbourke.net/dataformats/cube/)
  * Bader Charge Analysis: [http://theory.cm.utexas.edu/henkelman/code/bader](http://theory.cm.utexas.edu/henkelman/code/bader)


How it works:
    This program works in several steps:
    1. get middle single molecule of the constructed supercell, we will select
    hole positions according to this middle single molecule
    2. after you get the single molecule structure, please calculate the HOMO
    and put the .cube file inside /data/singlemolecule
    3. according to this singlemolecule HOMO, we will decide the Cartesian coordinates of holes and output the exciton wavefunction input file for BerkeleyGW
    4. after you run all the exciton wavefunction calculations, we will collect
    all charge files and run dba over and output the final charge transfer character

Functionality:  
    - find hole positions with HOMO cube file  
    - generate input file, which is define excition wfn calculation input file  
    - automate bader analysis over exciton wfn result and collect CT percentage  
    - check the convergence of exciton wfn calculation based on the output cube file and bader analysis

Reminder:
    - DO NOT rerun single molecule selection after you get the HOMO,
    because the order of this single molecule might change
    - Do rerun single molecule selection if you decide to change the
    supercell size

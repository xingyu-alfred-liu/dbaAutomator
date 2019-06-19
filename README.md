# dba_automation

If you are in a hurry, you can only read the contents within
"!!! IMPORTANT !!!"

!!! IMPORTANT !!!

Required Libraries:
    - ase
    - pymatgen
    - numpy

What you need to do:
    (1) put the mean-field calculation input file under the path: 'dba_automation/data/unitcell'
    (2) define the supercell size in the script: 'dba_automation/src/files.py' (recommended)
    OR
    put the fine grid file 'kgrid.in' into the path:
    'dba_automation/data/unitcell' (not recommended)
    (3) go to 'dba_automation/src' and run the script:
    'step1_selectSingleMolecule.py'

!!!!!!!!!!!!!!!!!

Introduction:
dba_automation is a simple method to automate double-bader analysis (DBA). dba_automation interfaces with the output of
excition wfn calculations (.cube files) performed by BerkeleyGW
(https://berkeleygw.org/), and most of the DFT codes cope with
Berkeley, such as Quantum ESPRESSO. 

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

Some external resources:  
    - cube file explanation: http://paulbourke.net/dataformats/cube/  
    - Bader implementation we used: http://theory.cm.utexas.edu/henkelman/code/bader
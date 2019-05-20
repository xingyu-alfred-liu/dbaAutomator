# dba_automation
A simple method to automate DBA process.

Functionality:  
    - find hole positions with HOMO cube file  
    - generate input file, which is define excition wfn calculation input file  
    - automate bader analysis over exciton wfn result and collect CT percentage  

Requirements:  
    - bader bin path  
    - HOMO charge path  
    - BerkeleyGW bin and exciton wfn path  

Some external resources:  
    - cube file explanation: http://paulbourke.net/dataformats/cube/  
    - Bader implementation we used: http://theory.cm.utexas.edu/henkelman/code/bader/  

DO NOT FORGET:  
    - cube file has Bohr and Angstrom units, FHI-aims use Bohr (at least for my output file), remember to convert it to Angstrom

"""
This is the file defining the path for necessary files,
please make sure the vars below are correct:
    
    - dataDir
    string
    which is the folder "data" included under this repo
    Inside data folder, there should be three subfolders:
    1. unitcell
    2. singlemolecule
    3. supercell
    Inside the unitcell, there should be an mean-field input 
    file, our default is using Quantum ESPRESSO but ase should 
    be able to handle most of the 

    - fineGrid
    list
    Define the find grid number, default will be an empty list
    unless user type in the fine grid setting themself
    Please make sure the first line of kgrid.in is the fine grid

    - baderPath
    string
    The path to bader analysis binary
    The Bader program can be downloaded here: 
    http://theory.cm.utexas.edu/henkelman/code/bader/
    By default, the bader analysis results will be output in the 
    binary path, so dba_automator will copy the bin and run it 
    in a separate path, after the calculation is done, the bin will
    be deleted.

    - chargePath
    string
    The path to single molecule HOMO cube file

    - singleMolPath
    string
    The path where the single molecule charge and bader analysis locates

"""

# for single molecule processing
fineGrid = []
dataDir = '/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data'
# for bader analysis calculation
baderPath = '/Users/alfredliu/Documents/Research/Software/Bader/macos/bader'
chargePath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/singleMolecule/eigen_density185.cube"
singleMolPath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/singleMolecule"
supercellPath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/supercell"
chargeThreshold = 0.01
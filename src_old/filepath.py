"""
This is the file defining the path or necessary settings,
please type in:
    
    - QEinputPath
    string
    The file path towards the mean field calculation input file
    The purpose for this variable is to provide the program the 
    unit cell so that we can construct a supercell
    The default format is Quantum ESPRESSO but any structure file
    compatible with ase.io should be file
    
    - fineGridpath
    string
    The fine grid file. kgrid.in
    The purpose for this variable is to find out the find grid setting
    in your calculation so that we can construct the exciton wfn 
    calculation input file. 
    Also, directly input it also works and is recommended. 
    
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
QEinputPath = '/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/crystal/in'
fineGridpath = '/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/crystal/kgrid.in'
fineGrid = []
dataDir = '/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data'
# for bader analysis calculation
baderPath = '/Users/alfredliu/Documents/Research/Software/Bader/macos/bader'
chargePath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/singleMolecule/eigen_density185.cube"
singleMolPath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/singleMolecule"
supercellPath = "/Users/alfredliu/Desktop/papersWorkingOn/dba_automation/data/supercell"
chargeThreshold = 0.01
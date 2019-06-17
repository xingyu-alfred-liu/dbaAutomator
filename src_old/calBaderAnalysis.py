import os
from filepath import *

def run_dba(singleMolPath, baderPath, chargeFile):
    os.system(baderPath+' '+chargeFile)
    print(os.listdir('./'))

if __name__ == "__main__":
    run_dba(singleMolPath, baderPath, chargePath)

import ase
import os

def get_geo(filepath):
    with open(, 'r') as inputfile, open('tmp.xyz', 'w') as outputfile:
    for line in inputfile:
        eachline = line.split()
        if eachline[0] == "6":
            element = "C"
        elif eachline[0] == "1":
            element = "H"
        outputfile.write(element+'   ')
        for i in range(len(eachline)-1):
            outputfile.write(str(float(eachline[i+1])*0.529177)+'   ')
        outputfile.write('\n') 


if __name__ == '__main__':
    path = '/home/xingyu/software/dba_automation/data/singleMol'
    filelist = os.listdir(path)
    for name in filelist:
        if name.endswith('.cube'):
            filepath = path+'/'+name
            break
    get_geo(filepath)

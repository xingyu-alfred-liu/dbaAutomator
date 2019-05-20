
with open('structure_out_HOMO.xyz', 'r') as inputfile, open('tmp.xyz', 'w') as outputfile:
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


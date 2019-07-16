from dbaAutomator.core import checker

path = ''
finegrid = []

dba = checker(path, finegrid)
dba.calprep()
dba.checkconv()

# Bader Analysis binary will be copied into this folder. 

We use the Bader analysis code from Henkelman's Group:
    http://theory.cm.utexas.edu/henkelman/code/bader/

Slight modifications are made in bader analysis:

    /bader/bader_mod.f90, line 1461:
    WRITE(100,'(1I5,4F12.6,2F13.6)') i,ions%r_car(i,:),bdr%ionchg(i),bdr%minsurfdist(i),bdr%ionvol(i)
    --->
    WRITE(100,'(1I5,4F20.6,2F13.6)') i,ions%r_car(i,:),bdr%ionchg(i),bdr%minsurfdist(i),bdr%ionvol(i)

such change is made becuase the output file "ACF.dat" sometimes could not output all numbers, which could lead 
to errors in loading the results from bader analysis. 
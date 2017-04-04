#!/bin/tcsh
#qsub -l low -cwd -hard -q "low.q@minos21" -pe mpi-fu 12 run.csh
echo "start program on"
mpirun $PWD/bin/solver

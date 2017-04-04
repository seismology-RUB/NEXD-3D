#!/bin/tcsh
#qsub -l low -l arch=lx24-amd64 -R Y -cwd -hard -q "low.q@minos17" -l os=Debian_5\*_lenny -pe mpi-fu 1 run.csh
#qsub -l low -cwd -hard -q "low.q@minos20" -pe mpi-fu 16 run.csh
setenv OMP_NUM_THREADS $NSLOTS
echo "start program on"
echo $NSLOTS
mpirun -np $NSLOTS $PWD/bin/solver

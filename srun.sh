#!/bin/bash
if make debug ; then
  nprocs=${1:-4}
  if [ "$nprocs" == "1" ] ; then
    echo ">>>Using 1 working unit"
  else
    echo ">>>Using $nprocs working units"
  fi
  mpirun -np $nprocs ./sdd_simple_parallel 2>&1 | tee test.dat
else
  make debug 2>&1 | tee lib.dat
  echo "*************Error at compilation ***************"
fi

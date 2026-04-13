#!/bin/bash

########First call (Petsc to nonPetsc)

ok_ns=0
ns=1
while read line 
do
   if (( ok_ns == 1 )); then
      ns=$line
      ok_ns=0
   fi
   if [ "$line" == "===Number of processors in meridian section (Input)" ]; then
      ok_ns=1
   fi
done < data_interpol
nproc=$(( ns ))

$2 $3$nproc $4 "./$1" petsc_nonpetsc  >> lis 2>&1



########Second call (NonPetsc to Petsc)

ok_ns=0
ns=1
while read line 
do
   if (( ok_ns == 1 )); then
      ns=$line
      ok_ns=0
   fi
   if [ "$line" == "===Number of processors in meridian section (Output)" ]; then
      ok_ns=1
   fi
done < data_interpol
nproc=$(( ns ))

$2 $3$nproc $4 "./$1" nonpetsc_petsc  >> lis 2>&1

#!/bin/bash
cp data_20 data
datatest=data

ok_ns=0
ok_nf=0
ns=1
nf=1
while read line
do
   if (( ok_ns == 1 )); then
      ns=$line
      ok_ns=0
   fi
   if (( ok_nf == 1 )); then
      nf=$line
      ok_nf=0
   fi
   if [ "$line" == "===Number of processors in meridian section" ]; then
      ok_ns=1
   fi
   if [ "$line" == "===Number of processors in Fourier space" ]; then
      ok_nf=1
   fi
done < "$datatest"
nproc=$(( ns*nf ))

cp suite_ns_S000_I001.Mesh_10_form.FEM suite_ns_S000_I002.Mesh_10_form.FEM
./mesh_interpolation.sh "../EXECUTABLE/a_mesh_interpolation_19_20.exe" "$1" "$2" "$4"
for FILE in suite_*; do mv "$FILE" "$(echo "$FILE" | sed 's/_I002//')"; done
cp regression_reference_20 regression_reference

$1 $2$nproc $4 ../EXECUTABLE/$3 regression
echo $?
#cp current_regression_reference current_regression_reference_20

#Clean up
rm -f fort.* Mesh_*_FE_* *.plt mesh_part* suite_* lis current_regression_reference data regression_reference

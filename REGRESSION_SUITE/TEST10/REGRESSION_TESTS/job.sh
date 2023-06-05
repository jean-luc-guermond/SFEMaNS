#!/bin/bash
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
     

$1 $2$nproc $4  ../EXECUTABLE/$3 regression
echo $?
#cp current_regression_reference regression_reference_10

#Clean up
rm -f fort.* Mesh_*_FE_* dd.plt mesh_part* current_regression_reference

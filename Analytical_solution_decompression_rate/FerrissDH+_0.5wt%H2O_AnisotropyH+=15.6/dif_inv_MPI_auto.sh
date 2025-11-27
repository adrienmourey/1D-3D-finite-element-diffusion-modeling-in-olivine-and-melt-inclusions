#!/bin/bash

Profiles="SW_7_2_1 SW_7_3_1 SW_7_3_2 SW_7_4_1 SW_7_6_1 SW_7_8_1 SW_7_10_1 K97_23_1_1_1 K97_23_4_1 K97_23_5_1 K97_23_8_1 K97_23_10_1 K97_23_16_1 K97_23_16_2 K97_23_16_3 K97_23_17_1" 


#SW_7_7_1

#Profiles="K97_23_15_1" 

for x in ${Profiles} ; do # Loop through relevant profiles

echo "${x}"

# CHAINS_4
f_dat="Lerner2024MIs_golden_pumice.csv" # Data file
f_degas="magma_sat_hawaii_0.5_short.csv" # degassing path 
f_out="${x}_MI_dpdt" # Name of output folder
#mv "$f_out"/* ~/.trash

mkdir $f_out
 
mpiexec -n 4 python3 H2O_mi_inversion_analytical.py $f_dat $f_degas ${x} $f_out"/vars_1" # Run on 4 cores

python3 multinest_marginals.py $f_out/vars_1 # Final plotting script

done





#!/bin/bash


Profiles="16IPE1BINT_Ol-17 16IPE1BINT_Ol-9 16IPE1BINT_Ol-1 16IPE1BINT_Ol-21 16IPE1BINT_Ol-28 16IPE1BINT_Ol-2 16IPE1BINT_Ol-39 16IPE1BINT_Ol-44 16IPE1BINT_Ol-47_Mar23 16IPE1BINT_Ol-48_2017 16IPE1BINT_Ol-48_Mar23  16IPE1BINT_Ol-58_2017 16IPE1BINT_Ol-58_Mar23 16IPE1BINT_Ol-5"

# 16IPE1BINT_Ol-57_Mar23 16IPE1BINT_Ol-56 unstable

for x in ${Profiles} ; do # Loop through relevant profiles

# copy files into the corresponding crystal folder

#cp ./KC_fO2.py ./${x}_pyMN/KC_fO2.py # fO2 module

cp ./DFENS_Ol_FeMg_1D_mT2_f.py ./${x}_pyMN/DFENS_Ol_FeMg_1D_mT2_f.py # DFENS script

cp ./multinest_marginals.py ./${x}_pyMN/multinest_marginals.py # plotting script / output information

cp ./pmc.py ./${x}_pyMN/pmc.py # PyMultiNest module files

cd ./${x}_pyMN # change directory to the crystal folder

echo "${x}"

# CHAINS_4
f_dat="${x}_XFo_filt_err_bse.csv" # Data file
f_modpar="${x}_modpar_lu.csv" # Model parameters for the crystal
f_angles="${x}_angles.csv" # Angles 
f_out='chains_FeMg_lu_MPI_2mT' # Name of output folder
#mv "$f_out"/* ~/.trash

mkdir $f_out
 
mpiexec -n 18 python3 DFENS_Ol_FeMg_1D_mT2_f.py $f_dat $f_modpar $f_angles $f_out"/vars_1" # Run DFENS on 20 cores using mpiexec and python3

python3 multinest_marginals.py $f_out/vars_1 # Final plotting script

cd ..

done





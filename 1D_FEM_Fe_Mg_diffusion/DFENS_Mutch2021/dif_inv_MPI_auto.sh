#!/bin/bash

Profiles="16IPE1BINT_Ol-10 16IPE1BINT_Ol-51 16IPE1BINT_Ol-13 16IPE1BINT_Ol-52 16IPE1BINT_Ol-16  16IPE1BINT_Ol-24 16IPE1BINT_Ol-7 16IPE1BINT_Ol-26 16IPE1BINT_Ol-8 16IPE1BINT_Ol-32 16IPE1BINT_Ol-46  16IPE1BINT_Ol-50"


#16IPE1BINT_Ol-4 16IPE1BINT_Ol-6a # Missing angles


for x in ${Profiles} ; do # Loop through relevant profiles

# copy files into the corresponding crystal folder

#cp ./KC_fO2.py ./${x}_pyMN/KC_fO2.py # fO2 module

cp ./DFENS_Ol_FeMg_1D_f.py ./${x}_pyMN/DFENS_Ol_FeMg_1D_f.py # DFENS script

cp ./multinest_marginals.py ./${x}_pyMN/multinest_marginals.py # plotting script / output information

cp ./pmc.py ./${x}_pyMN/pmc.py # PyMultiNest module files

cd ./${x}_pyMN # change directory to the crystal folder

echo "${x}"

# CHAINS_4
f_dat="${x}_XFo_filt_err_bse.csv" # Data file
f_modpar="${x}_modpar_lu.csv" # Model parameters for the crystal
f_angles="${x}_angles.csv" # Angles 
f_out='chains_FeMg_lu_MPI_fx' # Name of output folder
#mv "$f_out"/* ~/.trash

mkdir $f_out
 
mpiexec -n 17 python3 DFENS_Ol_FeMg_1D_f.py $f_dat $f_modpar $f_angles $f_out"/vars_1" # Run DFENS on 20 cores using mpiexec and python3

python3 multinest_marginals.py $f_out/vars_1 # Final plotting script

cd ..

done





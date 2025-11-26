#!/bin/bash

Profiles_1="Ol2B Ol13A"

Profiles_2="Ol5B_main Ol4A"

#16IPE1BINT_Ol-4 16IPE1BINT_Ol-6a # Missing angles


for x in ${Profiles_1} ; do # Loop through relevant profiles

# copy files into the corresponding crystal folder

#cp ./KC_fO2.py ./${x}_pyMN/KC_fO2.py # fO2 module

cp ./DFENS_Ol_FeMg_1D_f_plot.py ./${x}_pyMN/DFENS_Ol_FeMg_1D_f_plot.py # DFENS script

cp ./multinest_marginals.py ./${x}_pyMN/multinest_marginals.py # plotting script / output information

cp ./pmc.py ./${x}_pyMN/pmc.py # PyMultiNest module files

cd ./${x}_pyMN # change directory to the crystal folder

echo "${x}"

# CHAINS_4
f_dat="${x}_XFo_filt_err.csv" # Data file
f_modpar="${x}_modpar_lu.csv" # Model parameters for the crystal
f_angles="${x}_angles.csv" # Angles 
f_out='chains_FeMg_lu_MPI_fx' # Name of output folder
#mv "$f_out"/* ~/.trash

mkdir $f_out
 
mpiexec -n 5 python3 DFENS_Ol_FeMg_1D_f_plot.py $f_dat $f_modpar $f_angles $f_out"/vars_1" # Run DFENS on 20 cores using mpiexec and python3

python3 multinest_marginals.py $f_out/vars_1 # Final plotting script

cd ..

done



for x in ${Profiles_2} ; do # Loop through relevant profiles

# copy files into the corresponding crystal folder

#cp ./KC_fO2.py ./${x}_pyMN/KC_fO2.py # fO2 module

cp ./DFENS_Ol_FeMg_1D_mT2_plot_f.py ./${x}_pyMN/DFENS_Ol_FeMg_1D_mT2_plot_f.py # DFENS script

cp ./multinest_marginals.py ./${x}_pyMN/multinest_marginals.py # plotting script / output information

cp ./pmc.py ./${x}_pyMN/pmc.py # PyMultiNest module files

cd ./${x}_pyMN # change directory to the crystal folder

echo "${x}"

# CHAINS_4
f_dat="${x}_XFo_filt_err.csv" # Data file
f_modpar="${x}_modpar_lu.csv" # Model parameters for the crystal
f_angles="${x}_angles.csv" # Angles 
f_out='chains_FeMg_lu_MPI_2mT' # Name of output folder
#mv "$f_out"/* ~/.trash

mkdir $f_out
 
mpiexec -n 5 python3 DFENS_Ol_FeMg_1D_mT2_plot_f.py $f_dat $f_modpar $f_angles $f_out"/vars_1" # Run DFENS on 20 cores using mpiexec and python3

python3 multinest_marginals.py $f_out/vars_1 # Final plotting script

cd ..

done



#/bin/bash

# wrapper script

file_name1="olivine_2B" # Name of mesh files

temp=1210.0 #Temperature (C)
fo2=0.000887 #Oxygen fugacity in Pa
P_Pa=74000000.0 #Pressure in Pa
fo_core=0.844 #XFo content in olivine core
fo_rim=0.851 #XFo content in olivine rim
femg_core=0.426 #Fe/Mg ratio in olivine core
femg_rim=0.403 #Fe/Mg ratio in olivine rim
kd=0.335 #Olivine-melt partition coefficient for Fe-Mg
time=3000 #Days
nts=300 # Number of timesteps

mpirun -np 18 python3 3D_diffusion_FeMgm_ct_linear2.py ${file_name1} ${temp} ${fo2} ${P_Pa} ${fo_core} ${fo_rim} ${femg_core} ${femg_rim} ${kd} ${time} ${nts} # Currently being run on 18 cores



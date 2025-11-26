#/bin/bash

# wrapper script

file_name1="olivine_5B" # Name of mesh files

temp=1230.0 #diffusion anisotropy
fo2=0.0014922
P_Pa=74000000.0
fo_core=0.895
fo_rim=0.866
femg_core=0.269
femg_rim=0.355
kd=0.335 #Olivine-melt partition coefficient for water
time=11000
nts=550 # Number of timesteps

mpirun -np 18 python3 3D_diffusion_FeMgm_ct_linear2.py ${file_name1} ${temp} ${fo2} ${P_Pa} ${fo_core} ${fo_rim} ${femg_core} ${femg_rim} ${kd} ${time} ${nts} # Currently being run on 18 cores



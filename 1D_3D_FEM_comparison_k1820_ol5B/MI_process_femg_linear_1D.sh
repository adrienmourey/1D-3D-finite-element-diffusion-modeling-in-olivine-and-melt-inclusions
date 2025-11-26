#/bin/bash

# wrapper script

f_dat="Ol5B_3D_results.csv"
f_angle="Angles.csv"

temp=1230.0 #Temperature C
fo2=0.0014922 #fO2 Pa
P_Pa=74000000.0 #Pressure Pa
fo_core=0.895 #Forsterite composition core
fo_rim=0.866 #Forsterite composition rim
femg_core=0.269 # Fe/Mg core
femg_rim=0.355 # Fe/Mg rim
time=11000 # time days
nts=550 # Number of timesteps

python3 1D_diffusion_FeMgm_linear_comparison.py ${f_dat} ${f_angle} ${temp} ${fo2} ${P_Pa} ${fo_core} ${fo_rim} ${femg_core} ${femg_rim} ${time} ${nts} 


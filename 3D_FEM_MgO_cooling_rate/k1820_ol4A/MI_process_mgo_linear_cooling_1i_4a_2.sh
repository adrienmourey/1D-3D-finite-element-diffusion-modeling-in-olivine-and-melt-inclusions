#/bin/bash


files="olivine_4A"  #olivine_2B " # Name of mesh files

for x in ${files} ; do

echo ${x}

Ts=1213.0 # temperature at start
Te=1199.0 # temperature at end
fo2=0.000935
P_Pa=74000000.0
fo_av=0.875
kd=4.5 #Olivine-melt partition coefficient for Mg
nts=500 # Number of timesteps

cooling="0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5" # Cooling rate

for y in ${cooling} ; do

echo ${y}

cp ./3D_diffusion_MgO_ct_linear_cooling_1i.py ./${x}/3D_diffusion_MgO_ct_linear_cooling_1i.py

cd ./${x}

mpirun -np 17 python3 3D_diffusion_MgO_ct_linear_cooling_1i.py ${x} ${Ts} ${Te} ${fo2} ${P_Pa} ${fo_av} ${kd} ${y} ${nts} # Currently being run on 18 cores

cd ..

done

done

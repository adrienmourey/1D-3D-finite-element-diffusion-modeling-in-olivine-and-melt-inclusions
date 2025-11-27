#/bin/bash

# wrapper script

files="olivine_4A olivine_2B olivine_1A" # Name of mesh files
file_degas="Kilauea_closed_degassing_1210_vshort.csv"
file_melt="Lerner_melt.csv"


for x in ${files} ; do

echo ${x}

aniso=15.6 #diffusion anisotropy
nts=150 # Number of timesteps
kd=0.0009 #Olivine-melt partition coefficient for water
dpdt=0.45 # Magma decompression rate MPa/s
temp=1210.0 #Temperature (C)

echo ${dpdt}

cp ./3D_diffusion_Kilauea_ct_barth.py ./${x}/3D_diffusion_Kilauea_ct_barth.py
cp ./water_melt.py ./${x}/water_melt.py
cp ./${file_melt} ./${x}/${file_melt}
cp ./${file_degas} ./${x}/${file_degas}

cd ./${x}

mpirun -np 18 python3 3D_diffusion_Kilauea_ct_barth.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

cd ..

done



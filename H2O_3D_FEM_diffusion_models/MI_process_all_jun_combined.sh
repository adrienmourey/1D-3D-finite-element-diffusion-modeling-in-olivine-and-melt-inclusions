#/bin/bash

# wrapper script

#files=" olivine_1A" # Name of mesh files
files="olivine_4A olivine_2B" # Name of mesh files olivine_4A olivine_2B
file_degas="magma_sat_hawaii_0.6_short.csv"
file_melt="Lerner_melt.csv"

files_2="olivine_1A"

barth="Barth"
ferris="Ferris"

for x in ${files} ; do

echo ${x}

aniso=15.6 #diffusion anisotropy
nts=10000 # Number of timesteps
kd=0.0009 #Olivine-melt partition coefficient for water
dpdt=0.07 # Magma decompression rate MPa/s
temp=1210.0 #Temperature (C)

echo ${dpdt}

cp 3D_diffusion_Kilauea_ct_barth_et.py ./${x}/3D_diffusion_Kilauea_ct_barth_et.py
cp 3D_diffusion_Kilauea_ct_ferris_et.py ./${x}/3D_diffusion_Kilauea_ct_ferris_et.py
cp water_melt.py ./${x}/water_melt.py
cp ${file_melt} ./${x}/${file_melt}
cp ${file_degas} ./${x}/${file_degas}

cd ./${x}

echo ${barth}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_barth_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

echo ${ferris}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_ferris_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores


aniso=15.6 #diffusion anisotropy
nts=10000 # Number of timesteps
kd=0.0009 #Olivine-melt partition coefficient for water
dpdt=0.001 # Magma decompression rate MPa/s
temp=1210.0 #Temperature (C)

echo ${dpdt}

echo ${barth}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_barth_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

echo ${ferris}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_ferris_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

aniso=15.6 #diffusion anisotropy
nts=10000 # Number of timesteps
kd=0.0009 #Olivine-melt partition coefficient for water
dpdt=0.45 # Magma decompression rate MPa/s
temp=1210.0 #Temperature (C)

echo ${dpdt}

echo ${barth}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_barth_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

echo ${ferris}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_ferris_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

cd ..

done


for x in ${files_2} ; do

echo ${x}

aniso=15.6 #diffusion anisotropy
nts=15000 # Number of timesteps
kd=0.0009 #Olivine-melt partition coefficient for water
dpdt=0.45 # Magma decompression rate MPa/s
temp=1210.0 #Temperature (C)

echo ${dpdt}

cp 3D_diffusion_Kilauea_ct_barth_et.py ./${x}/3D_diffusion_Kilauea_ct_barth_et.py
cp 3D_diffusion_Kilauea_ct_ferris_et.py ./${x}/3D_diffusion_Kilauea_ct_ferris_et.py
cp water_melt.py ./${x}/water_melt.py
cp ${file_melt} ./${x}/${file_melt}
cp ${file_degas} ./${x}/${file_degas}

cd ./${x}

echo ${barth}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_barth_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

echo ${ferris}
mpirun -np 18 python3 3D_diffusion_Kilauea_ct_ferris_et.py ${x} ${file_melt} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${temp} # Currently being run on 18 cores

cd ..

done







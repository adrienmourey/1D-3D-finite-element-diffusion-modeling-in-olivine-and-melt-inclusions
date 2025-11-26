#!/bin/bash

# script to collate pymultinest outputs into one folder for processing. 

# Titan01 profiles
Profiles_1="Ol1A Ol2B Ol5B_lg Ol13A"

Profiles_2="Ol5B_main Ol4A"

rm Ol_1step_pyMN_time_stats.csv

touch Ol_1step_pyMN_time_stats.csv

rm Ol_1step_pyMN_medians.csv

touch Ol_1step_pyMN_medians.csv

rm Ol_2step_pyMN_time_stats.csv

touch Ol_2step_pyMN_time_stats.csv

rm Ol_2step_pyMN_medians.csv

touch Ol_2step_pyMN_medians.csv

for x in ${Profiles_1} ; do 

echo "${x}"

f_out='chains_FeMg_lu_MPI_fx' # Name of output folder

#cp ./multinest_marginals_blue.py ./${x}_pyMN_output/multinest_marginals_blue.py

cp ./Multinest_stats_out.py ./${x}_pyMN/Multinest_stats_out.py

cp ./Multinest_dfout.py ./${x}_pyMN/Multinest_dfout.py

cd ${x}_pyMN

#python multinest_marginals_blue.py vars_1 5

python3 Multinest_dfout.py $f_out/vars_1 ${x} # Final plotting script

python3 Multinest_stats_out.py $f_out/vars_1  ${x}

cd ..

cat ./${x}_pyMN/${x}_pyMN_stats.csv >> ./Ol_1step_pyMN_time_stats.csv

cat ./${x}_pyMN/${x}_pyMN_medians.csv >> ./Ol_1step_pyMN_medians.csv

done


for x in ${Profiles_2} ; do 

echo "${x}"

f_out='chains_FeMg_lu_MPI_2mT' # Name of output folder

#cp ./multinest_marginals_blue.py ./${x}_pyMN_output/multinest_marginals_blue.py

cp ./Multinest_stats_out_multi.py ./${x}_pyMN/Multinest_stats_out_multi.py

cp ./Multinest_dfout_multi.py ./${x}_pyMN/Multinest_dfout_multi.py

cd ${x}_pyMN

#python multinest_marginals_blue.py vars_1 5

python3 Multinest_dfout_multi.py $f_out/vars_1 ${x} # Final plotting script

python3 Multinest_stats_out_multi.py $f_out/vars_1  ${x}

cd ..

cat ./${x}_pyMN/${x}_pyMN_stats.csv >> ./Ol_2step_pyMN_time_stats.csv

cat ./${x}_pyMN/${x}_pyMN_medians.csv >> ./Ol_2step_pyMN_medians.csv

done

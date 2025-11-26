#!/bin/bash

# script to collate pymultinest outputs into one folder for processing. 

# Titan01 profiles
Profiles="LS_2156a1B_XENOLA_P1 LS_2156F_OLG_P2 LS_2156F_OLJ_P1 LS_2156F_XENOLE_P1 LS_2156a1B_2_XENOLC_P1 LS_2156a1B_OLA_P1 LS_2156a1B_OLA_P2 LS_2156a1B_2_OLD_P1 LS_2156a1B_2_XENOLB_P1 LS_2156F_OLG_P1"

#SKU_4_C2_P1 KRA84-10_cpxtr KRA9-1_cpxtr

rm LS_1step_pyMN_time_stats.csv

touch LS_1step_pyMN_time_stats.csv

rm LS_1step_pyMN_medians.csv

touch LS_1step_pyMN_medians.csv

for x in ${Profiles} ; do

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

cat ./${x}_pyMN/${x}_pyMN_stats.csv >> ./LS_1step_pyMN_time_stats.csv

cat ./${x}_pyMN/${x}_pyMN_medians.csv >> ./LS_1step_pyMN_medians.csv

done














































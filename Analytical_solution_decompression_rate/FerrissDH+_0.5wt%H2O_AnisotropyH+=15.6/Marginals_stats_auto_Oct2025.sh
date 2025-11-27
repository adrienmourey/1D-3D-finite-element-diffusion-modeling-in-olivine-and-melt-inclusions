#!/bin/bash

# script to collate pymultinest outputs into one folder for processing. 

# Titan01 profiles
Profiles="SW_7_2_1 SW_7_3_1 SW_7_3_2 SW_7_4_1 SW_7_6_1 SW_7_7_1  SW_7_10_1 K97_23_4_1 K97_23_5_1 K97_23_8_1 K97_23_10_1 K97_23_16_3 K97_23_17_1" 

# K97_23_16_1 K97_23_16_2 SW_7_8_1 K97_23_1_1_1

rm Kil_pyMN_time_stats.csv

touch Kil_pyMN_time_stats.csv

rm Kil_pyMN_medians.csv

touch Kil_pyMN_medians.csv

rm Kil_pyMN_dfout.csv

touch Kil_pyMN_dfout.csv

for x in ${Profiles} ; do

echo "${x}"

f_out="${x}_MI_dpdt"

#cp ./multinest_marginals_blue.py ./${x}_pyMN_output/multinest_marginals_blue.py

#python multinest_marginals_blue.py vars_1 5

python3 Multinest_dfout.py $f_out/vars_1 ${x} # Final plotting script

python3 Multinest_stats_out.py $f_out/vars_1  ${x}

cat ./${x}_pyMN_stats.csv >> ./Kil_pyMN_time_stats.csv

cat ./${x}_pyMN_medians.csv >> ./Kil_pyMN_medians.csv

cat ./${x}_pyMN_dataframe.csv >> ./Kil_pyMN_dfout.csv


mv ./${x}_pyMN_stats.csv ./$f_out/${x}_pyMN_stats.csv

mv ./${x}_pyMN_medians.csv ./$f_out/${x}_pyMN_medians.csv

mv ./${x}_pyMN_dataframe.csv ./$f_out/${x}_pyMN_dataframe.csv

mv ./${x}_CFD.pdf ./$f_out/${x}_CFD.pdf


# Combine df_out as well


done














































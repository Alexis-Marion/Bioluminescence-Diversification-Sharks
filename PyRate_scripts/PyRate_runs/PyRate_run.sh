# This script run in parallel several pyrate runs from the same set of replicates
# Assuming you have parallel installed 20 CPU will be used
# Occ_gn_Squali_PyRate.py is the location of your PyRate occurence file (in this case the occurence for genus data including singletons)
# ../../PyRate/PyRate.py is the location of your PyRate launch file
# epochs_rate_shift_5_Myrs.txt is your epoch file for shifts (truncated 5 Myrs in this case)
# Preservation/epochs.txt is your epoch file for preservation


echo "python ../../PyRate/PyRate.py  Data/Occ_gn_Squali_PyRate.py -fixShift epochs_rate_shift_5_Myrs.txt -j \$1 -qShift Preservation/epochs.txt -mG -n 20000000 -s 20000 -N 23" > tmp_script.sh
		parallel -j 20 bash tmp_script.sh ::: {1..20}

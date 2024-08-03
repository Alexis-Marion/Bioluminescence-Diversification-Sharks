# This is the script used to perform the automatic selection for the best-fit preservation model
# ../../PyRate/PyRate.py designate the location of your PyRate.py file
# $1 Is the location of your PyRate occurence file
# Assuming you have parallel installed, $2 is the number of CPU you want to allow for computation
# epochs.txt is the file with epochs that will be tested using PPmodeltest from PyRate
# model_drafting.r is the .r file used to perform model selection
# Other files are automatically generated

cntr=0
declare -i cntr
while [ ! -f result_boxplot.pdf ]; do
		cntr+=1
		echo "$cntr round of optimization"
		echo "python ../../PyRate/PyRate.py  $1 -j \$1 -qShift epochs.txt -mG -PPmodeltest > File\$1.log" > tmp_script.sh
		parallel -j $2 bash tmp_script.sh ::: {1..100}
		bash script.sh $cntr
		Rscript model_drafting.r tab.tsv epochs.txt
		rm tab.tsv
		rm tmp_script.sh
		sleep 1
done

# This is the master script for plotting the results from a PyRate run
# $1 where the output of a simple PyRate run
# $2 a txt file with the geological ages for the preservation rates
# $3 the number of burnin state you want to discard (see assess_run_convergence.py)
# $4 an output directory for the plotting
# $5 a txt file with the time bins
# $6 type of analysis (RJMCMC/BDS)
# All other arguments are either generated by this script or are associated with the folder for this analysis

mkdir $4/

## Ltt

mkdir $4/output_ltt
bash ltt_creator.sh $1 $3 $4/output_ltt

## Preservation rates

bash Q_rate.sh $1 $3 $4 $2

## Rtt & Plot

if [$6 == "BDS"]
then
	python ../../PyRate/PyRate.py -plot $1 -tag KEEP -b $3
	for file in $1/*rate_RTT.r
	do
		Rscript Plot_rates.r $file $5 $2 $4/output_ltt/ $4/Parsed_Q_rates.tsv $4/All_in_one_plot.pdf $6
	done
else
	python ../../PyRate/PyRate.py -plot $1 -tag KEEP -b $3
	for file in $1/*rates_RTT.r
	do
		Rscript Plot_rates.r $file $5 $2 $4/output_ltt/ $4/Parsed_Q_rates.tsv $4/All_in_one_plot.pdf $6
	done
fi
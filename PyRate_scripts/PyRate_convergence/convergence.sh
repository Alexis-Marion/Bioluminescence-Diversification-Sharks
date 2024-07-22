# assess_run_convergence.py is the python script computing ESS value for a set of pyrate log files
# ../Dummy_datapyrate_mcmc_logs/ is the path to the PyRate log file
# BDS is the type of analyses ran (Either RJMCMC or BDS). Here this value is BDS, namely birth-death with time-constrained shifts
# plot_ess.r is a plotting script used for representing the posterior probability ESS value and other useful metrics extracted from the predifined set of PyRate log files
# Each file with posterior proability above 200 is flagged as "KEEP", and will be selected for further analyses

python assess_run_convergence.py -dir ../Dummy_data/pyrate_mcmc_logs/ -ana BDS

Rscript plot_ess.r ../Dummy_data/pyrate_mcmc_logs/ESS_summary.txt Stat_PyRate_run.pdf
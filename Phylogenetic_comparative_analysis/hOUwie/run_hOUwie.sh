# This script runs multiple hOUWie analyses and assumes you have 20 CPU available
# Two scripts are used here, one for the consensus tree (hOUwie_consensus.r) and another one for the posterior tree distribution (hOUwie_replicated.r)
# The first argument argument is the trait examined (7 for bioluminescence; 9 for habitat)
# The second/last argument is the prefix for the output directory/file

# Consensus tree

# Bioluminescence

Rscript hOUwie_consensus.r 7 bioluminescence

# Habitat

Rscript hOUwie_consensus.r 9 habitat

# Posterior distribution

# Bioluminescence

echo "Rscript Panova_hOUwie_consensus_habitat_Replicated.r \$1 7 bioluminescence" > tmp_script.sh
		parallel -j 20 bash tmp_script.sh ::: {1..100}

# Habitat

echo "Rscript Panova_hOUwie_consensus_habitat_Replicated.r \$1 9 habitat" > tmp_script.sh
		parallel -j 20 bash tmp_script.sh ::: {1..100}
# This script runs multiple corHMM analyses
# Two scripts are used here, one for the consensus tree (OUwie_consensus.r) and another one for the posterior tree distribution (OUwie_replicated.r)
# The first argument is the path to the corHMM anaylis used for mapping the trait on the phylogeny
# The second argument is the trait examined (7 for bioluminescence; 9 for habitat)
# The third argument is the path to 
# The last argument is the prefix for the output directory/file

# Consensus tree

# Bioluminescence

Rscript OUwie_consensus.r ../corHMM/corHMM/bioluminescence.rds 7 OUwie/bioluminescence.tsv OUwie/bioluminescence.rds

# Bioluminescence with fossil

Rscript OUwie_consensus.r ../corHMM/corHMM/bioluminescence_fossil.rds 7 OUwie/bioluminescence_fossil.tsv OUwie/bioluminescence_fossil.rds

# Habitat

Rscript OUwie_consensus.r ../corHMM/corHMM/habitat.rds 9 OUwie/habitat.tsv OUwie/habitat.rds

# Habitat with fossil

Rscript OUwie_consensus.r ../corHMM/corHMM/habitat_fossil.rds 9 OUwie/habitat_fossil.tsv OUwie/habitat_fossil.rds

# Posterior distribution

# Bioluminescence

Rscript OUwie_replicated.r ../corHMM/corHMM/bioluminescence/ 7 OUwie/bioluminescence.tsv OUwie/bioluminescence/

# Bioluminescence with fossil

Rscript OUwie_replicated.r ../corHMM/corHMM/bioluminescence_fossil/ 7 OUwie/bioluminescence_fossil.tsv OUwie/bioluminescence_fossil/

# Habitat

Rscript OUwie_replicated.r ../corHMM/corHMM/habitat/ 9 OUwie/habitat.tsv OUwie/habitat/

# Habitat with fossil

Rscript OUwie_replicated.r ../corHMM/corHMM/habitat_fossil/ 9 OUwie/habitat_fossil.tsv OUwie/habitat_fossil/


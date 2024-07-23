# This script runs multiple corHMM analyses
# Two scripts are used here, one for the consensus tree (corHMM_ASE_consensus.r) and another one for the posterior tree distribution (corHMM_ASE_replicated.r)
# The first argument is the phylogeny used (Extant only/Including fossils)
# The second argument is the trait examined (3 for bioluminescence; 4 for habitat)
# The last argument is the prefix for the output directory/file

# Consensus tree

# Bioluminescence

Rscript corHMM_ASE_consensus.r ../../raw_data/Squaliformes_extant.tree 3 bioluminescence

# Bioluminescence with fossil

Rscript corHMM_ASE_consensus.r ../../raw_data/Squaliformes_fossil.tree 3 bioluminescence_fossil

# Habitat

Rscript corHMM_ASE_consensus.r ../../raw_data/Squaliformes_extant.tree 4 habitat

# Habitat with fossil

Rscript corHMM_ASE_consensus.r ../../raw_data/Squaliformes_fossil.tree 4 habitat_fossil

# Posterior distribution

# Bioluminescence

Rscript corHMM_ASE_replicated.r ../../raw_data/Squaliformes_posterior_distribution.tree 3 bioluminescence

# Bioluminescence with fossil

Rscript corHMM_ASE_replicated.r ../../raw_data/Multi_Squaliformes_fossil_posterior_distribution.tree 3 bioluminescence_fossil

# Habitat

Rscript corHMM_ASE_replicated.r ../../raw_data/Squaliformes_posterior_distribution.tree 4 habitat

# Habitat with fossil

Rscript corHMM_ASE_replicated.r ../../raw_data/Multi_Squaliformes_fossil_posterior_distribution.tree 4 habitat_fossil
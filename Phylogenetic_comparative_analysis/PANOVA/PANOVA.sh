# This scripts runs several Phylogenetic Analysis of Variance (PANOVA)
# Two subscripts will perform PANOVA on the consensus tree (Phylogenetic_analysis_of_variance_(PANOVA).r) and posterior distribution (Rscript Phylogenetic_analysis_of_variance_(PANOVA-Replicated).r)
# Two analysis for each dataset are performed, one for bioluminescence and an other for habitat

# Run the PANOVA on the consensus tree

# Bioluminescence

Rscript Phylogenetic_analysis_of_variance_(PANOVA).r  7 bioluminescence

# Habitat

Rscript Phylogenetic_analysis_of_variance_(PANOVA).r  9 habitat

# Run the PANOVA on the posterior tree distribution

# Bioluminescence

echo "Rscript Phylogenetic_analysis_of_variance_(PANOVA-Replicated).r \$1 7 bioluminescence" > tmp_biol_script.sh
		parallel -j 40 bash tmp_biol_script.sh ::: {1..100}

# Habitat

echo "Rscript Phylogenetic_analysis_of_variance_(PANOVA-Replicated).r \$1 9 habitat" > tmp_biol_script.sh
		parallel -j 40 bash tmp_biol_script.sh ::: {1..100}

library("stringr")
library("phytools")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.tree("../../raw_data/Squaliformes_extant.tree")

df$Species<-str_replace(df$Species, " ", "_")
states<-df
states_traits<-states[!states$Species %in% setdiff(states$Species, phy$tip.label),]

categorical_trait<-states_traits[args[1]]
names(categorical_trait)<-states_traits$Species
continuous_trait<-log(states_traits$Body.size)
names(continuous_trait)<-states_traits$Species

phanov<-phylANOVA(phy, categorical_trait, continuous_trait, nsim=10000, posthoc=TRUE, p.adj="bonferroni")

saveRDS(phanov, paste("PANOVA/data_panova_", args[2], "_", ".rds", sep = '_'))

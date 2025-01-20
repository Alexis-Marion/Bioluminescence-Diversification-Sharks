library("stringr")
library("phytools")
library("RRPP")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.tree("../../raw_data/Squaliformes_extant.tree")

df$Species<-str_replace(df$Species, " ", "_")
states<-df
states_traits<-states[!states$Species %in% setdiff(states$Species, phy$tip.label),]

categorical_trait<-states_traits[,as.numeric(args[1])]
names(categorical_trait)<-states_traits$Species
continuous_trait<-log(states_traits$Body.size)
names(continuous_trait)<-states_traits$Species

phanov_sim <- phylANOVA(phy, categorical_trait, continuous_trait, nsim=10000, posthoc=TRUE, p.adj="bonferroni")

GLS_PANOVA_RRPP <- lm.rrpp(continuous_trait ~ categorical_trait, 
                  print.progress = FALSE,
                  Cov = vcv(phy),
                  turbo = FALSE, verbose = TRUE)

saveRDS(phanov_sim, paste("panova_sim_", args[2], ".rds", sep = ""))

saveRDS(anova(GLS_PANOVA_RRPP), paste("panova_RRPP_", args[2], ".rds", sep = ""))

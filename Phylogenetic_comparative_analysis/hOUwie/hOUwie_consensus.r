library("stringr")
library("phytools")
library("ggstatsplot")
library("tidyverse")
library("OUwie")
library("geiger")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.tree("../../raw_data/Squaliformes_extant.tree")

df$Species<-str_replace(df$Species, " ", "_")
states<-df
states_traits<-states[!states$Species %in% setdiff(states$Species, phy$tip.label),
categorical_trait<-states_traits[,args[1]]
continuous_trait<-log(states_traits$Body.size)
data<-as.data.frame(cbind(states_traits$Species, categorical_trait,continuous_trait))
data$continuous_trait<-as.numeric(data$continuous_trait)

R1_bm1_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "BM1", rate.cat = 1, nSim = 100)

R1_ou1_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "OU1", rate.cat = 1, nSim = 100)

R1_oum_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "OUM", rate.cat = 1, nSim = 100)

R1_bm_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "BM1", rate.cat = 1, nSim = 100)

R1_ou_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "OU1", rate.cat = 1, nSim = 100)

R1_oum_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "OUM", rate.cat = 1, nSim = 100)

R2_bm_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "BM1", rate.cat = 2, nSim = 100)

R2_ou_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "OU1", rate.cat = 2, nSim = 100)

R2_oum_er <- hOUwie(phy, data, discrete_model = "ER", continuous_model = "OUM", rate.cat = 2, nSim = 100)

R2_bm_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "BM1", rate.cat = 2, nSim = 100)

R2_ou_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "OU1", rate.cat = 2, nSim = 100)

R2_oum_ard <- hOUwie(phy, data, discrete_model = "ARD", continuous_model = "OU1", rate.cat = 2, nSim = 100)

model_set <- list(R1_bm_er = R1_bm_er, R1_ou_er = R1_ou_er, R1_oum_er = R1_oum_er,
         R1_bm_ard = R1_bm_ard, R1_ou_ard = R1_ou_ard, R1_oum_ard = R1_oum_ard,
         R2_bm_er = R2_bm_er, R2_ou_er = R2_ou_er, R2_oum_er = R2_oum_er,
         R2_bm_ard = R2_bm_ard, R2_ou_ard = R2_ou_ard, R2_oum_ard = R2_oum_ard)
print(getModelTable(model_set, type = "AICc"))

write.table(getModelTable(model_set, type = "AICc"), paste("hOUwie/table_hOUwie_", args[2], ".tsv", sep = '_'), sep ="\t")

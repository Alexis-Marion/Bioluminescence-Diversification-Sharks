library("stringr")
library("phytools")
library("ggstatsplot")
library("tidyverse")
library("OUwie")
library("geiger")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.tree("../../raw_data/Squaliformes_posterior_distribution.tree")
phy<-phy[[as.numeric(args[1])]]

df$Species<-str_replace(df$Species, " ", "_")

states<-df

states_traits<-states[!states$Species %in% setdiff(states$Species, phy$tip.label),]

categorical_trait<-states_traits$Habitat
names(categorical_trait)<-states_traits$Species

continuous_trait<-log(states_traits$Body.size)
names(continuous_trait)<-states_traits$Species

dtf<-as.data.frame(cbind(categorical_trait,continuous_trait))
dtf$continuous_trait<-as.numeric(dtf$continuous_trait)

dat<-cbind(rownames(dtf),dtf)

dat$continuous_trait<-as.numeric(dat$continuous_trait)

CID_OUM_structure<-getOUParamStructure(model = "OUM", nObsState = 2, rate.cat = 2, null.model = FALSE)
CID_OUM_structure[3,2]<-3
CID_OUM_structure[3,c(3,4)]<-4
print(CID_OUM_structure)

R1_bm1_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "BM1", rate.cat = 1, nSim = 100)
R1_ou1_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "OU1", rate.cat = 1, nSim = 100)
R1_oum_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "OUM", rate.cat = 1, nSim = 100)

R1_bm1_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "BM1", rate.cat = 1, nSim = 100)
R1_ou1_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "OU1", rate.cat = 1, nSim = 100)
R1_oum_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "OUM", rate.cat = 1, nSim = 100)

R2_bm1_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "BM1", rate.cat = 2, nSim = 100)
R2_ou1_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "OU1", rate.cat = 2, nSim = 100)
R2_oum_er <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = "OUM", rate.cat = 2, nSim = 100)

R2_bm1_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "BM1", rate.cat = 2, nSim = 100)
R2_ou1_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "OU1", rate.cat = 2, nSim = 100)
R2_oum_ard <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = "OUM", rate.cat = 2, nSim = 100)

R2_oum_er_CID <- hOUwie(phy, dat, discrete_model = "ER", continuous_model = CID_OUM_structure, rate.cat = 2, nSim = 100)
R2_oum_ard_CID <- hOUwie(phy, dat, discrete_model = "ARD", continuous_model = CID_OUM_structure, rate.cat = 2, nSim = 100)

model_set <- list(R1_bm1_er = R1_bm1_er, R1_ou1_er = R1_ou1_er, R1_oum_er = R1_oum_er,
         R1_bm1_ard = R1_bm1_ard, R1_ou1_ard = R1_ou1_ard, R1_oum_ard = R1_oum_ard,
         R2_bm1_er = R2_bm1_er, R2_ou1_er = R2_ou1_er, R2_oum_er = R2_oum_er,
         R2_bm1_ard = R2_bm1_ard, R2_ou1_ard = R2_ou1_ard, R2_oum_ard = R2_oum_ard,
         R2_oum_er_CID = R2_oum_er_CID, R2_oum_ard_CID = R2_oum_ard_CID)
print(getModelTable(model_set, type = "AICc"))

write.table(getModelTable(model_set, type = "AICc"), paste("hOUwie/table_hOUwie_hab", args[1], ".tsv", sep = '_'), sep ="\t")

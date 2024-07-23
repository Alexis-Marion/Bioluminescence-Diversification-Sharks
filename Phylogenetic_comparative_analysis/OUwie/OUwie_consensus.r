library("stringr")
library("phytools")
library("ggstatsplot")
library("tidyverse")
library("OUwie")
library("geiger")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files

fossil_list<-c("Protosqualus_albertsi", "Protosqualus_sigei", "Squalus_minor", "Protocentrophorus_balticus", "Squaliodalatias_savoiei", "Eodalatias_crenulatus", "Isistius_trituratus", "Cretascymnus_adonis", "Trigonognathus_virginiae", "Proetmopterus_hemmooriensis", "Eoetmopterus_supracretaceus", "Cretascymnus_westfalicus")

    phy<-readRDS(args[1])$phy
    node_lab<- phy$node.label
    rate_cat <-readRDS(args[1])$rate.cat
    eq_table<-unique(readRDS(args[1])$data[,2])
    if (rate_cat == 1){
        node_lab<- gsub("1", eq_table[1], node_lab)
        node_lab<- gsub("2", eq_table[2], node_lab)
    }
    if (rate_cat == 2){
        node_lab<- gsub("1", eq_table[1], node_lab)
        node_lab<- gsub("2", eq_table[2], node_lab)
        node_lab<- gsub("3", eq_table[1], node_lab)
        node_lab<- gsub("4", eq_table[2], node_lab)
    }
    if (rate_cat == 3){
        node_lab<- gsub("1", eq_table[1], node_lab)
        node_lab<- gsub("2", eq_table[2], node_lab)
        node_lab<- gsub("3", eq_table[1], node_lab)
        node_lab<- gsub("4", eq_table[2], node_lab)
        node_lab<- gsub("5", eq_table[1], node_lab)
        node_lab<- gsub("6", eq_table[2], node_lab)
    }
    phy$node.label<-node_lab
    phy<-drop.tip(phy, fossil_list)

df$Species<-str_replace(df$Species, " ", "_")

states<-df

states_traits<-states[!states$Species %in% setdiff(states$Species, phy$tip.label),]

categorical_trait<-states_traits[,as.numeric(args[2])]
names(categorical_trait)<-states_traits$Species

continuous_trait<-log(states_traits[,6])
names(continuous_trait)<-states_traits$Species

dtf<-as.data.frame(cbind(categorical_trait,continuous_trait))

dtf$continuous_trait<-as.numeric(dtf$continuous_trait)

dat<-cbind(rownames(dtf),dtf)

BM1 <- OUwie(phy, dat, model="BM1", algorithm = 'invert', simmap.tree = FALSE)
OU1 <- OUwie(phy, dat, model="OU1", algorithm = 'invert', simmap.tree = FALSE)
OUM <- OUwie(phy, dat, model="OUM", algorithm = 'invert', simmap.tree = FALSE)

np<-c(BM1$param.count, OU1$param.count, OUM$param.count)
lnLik<-c(BM1$loglik, OU1$loglik, OUM$loglik)
AICc<-c(BM1$AICc, OU1$AICc, OUM$AICc)
dAICc<-geiger::aicw(AICc)[,2]
AICcwt<-geiger::aicw(AICc)[,3]

model_table_AICc<-as.data.frame(cbind(np,lnLik, AICc, dAICc, AICcwt))
rownames(model_table_AICc)<-c("BM1", "OU1", "OUM")

write.table(model_table_AICc, args[3], sep ="\t")
saveRDS(eval(parse(text = rownames(model_table_AICc)[which.min(model_table_AICc$AICc)])),args[4])

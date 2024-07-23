library("corHMM")
library("stringr")
library("phytools")
library("qpcR")

args<-commandArgs(trailingOnly = TRUE)

df<-read.csv("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep="\t") # omit sep ="\t" for .csv files
phy<-read.tree(args[1])

df$Species<-str_replace(df$Species, " ", "_")

states<-cbind(df$Species, df$Body.size, df$Bioluminescent, df$Habitat)

states_traits<-states[!states[,1] %in% setdiff(states[,1], phy$tip.label),]

states_traits[is.na(states_traits)]<-"?"

states_traits<-states_traits[,c(1, as.numeric(args[2]))]
colnames(states_traits)<-c("Species", "Cat")

trait_1_eq<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

### Model 2 : all rates differ

trait_1_ard<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

### Model 1 : equal rates

trait_2_eq<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

### Model 2 : all rates differ

trait_2_ard<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

### Model 1 : equal rates

trait_3_eq<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

### Model 2 : all rates differ

trait_3_ard<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = TRUE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

data_trait<-data.frame(
                  cbind(c(trait_1_eq$loglik, trait_1_ard$loglik,
                          trait_2_eq$loglik, trait_2_ard$loglik,
                          trait_3_eq$loglik, trait_3_ard$loglik),
                        c(trait_1_eq$AICc, trait_1_ard$AICc,
                          trait_2_eq$AICc, trait_2_ard$AICc,
                          trait_3_eq$AICc, trait_3_ard$AICc),
                akaike.weights(c(trait_1_eq$AICc, trait_1_ard$AICc,
                          trait_2_eq$AICc, trait_2_ard$AICc,
                          trait_3_eq$AICc,  trait_3_ard$AICc))$deltaAIC,
                akaike.weights(c(trait_1_eq$AICc, trait_1_ard$AICc,
                          trait_2_eq$AICc, trait_2_ard$AICc, 
                          trait_3_eq$AICc, trait_3_ard$AICc))$weights,
                c((max(as.vector(trait_1_eq$index.mat)[!is.na(as.vector(trait_1_eq$index.mat))])), (max(as.vector(trait_1_ard$index.mat)[!is.na(as.vector(trait_1_ard$index.mat))])),
(max(as.vector(trait_2_eq$index.mat)[!is.na(as.vector(trait_2_eq$index.mat))])), (max(as.vector(trait_2_ard$index.mat)[!is.na(as.vector(trait_2_ard$index.mat))])),
(max(as.vector(trait_3_eq$index.mat)[!is.na(as.vector(trait_3_eq$index.mat))])), (max(as.vector(trait_3_ard$index.mat)[!is.na(as.vector(trait_3_ard$index.mat))])))
                ))
rownames(data_trait)<-c("trait_1_eq", "trait_1_ard",
                        "trait_2_eq",  "trait_2_ard",
                        "trait_3_eq",  "trait_3_ard")
colnames(data_trait)<-c("Log-likelihood", "AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(data_trait, paste("corHMM/", args[3], ".tsv", sep =""), sep ="\t") 

saveRDS(eval(parse(text = rownames(data_trait)[which.min(data_trait$AICc)])),paste("corHMM/", args[3], ".rds", sep =""))

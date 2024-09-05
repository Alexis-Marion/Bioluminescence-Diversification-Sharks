library("ape")

trees<-read.nexus("path_to_the posterior_tree_distribution")

set.seed(123)
phylo_list<-list()
class(phylo_list)<-"multiPhylo"
for (i in 1:100){
    sample_data<-sample(10000, 100)
    temp_tree<-trees[[sample_data[1]]]
    temp_tree<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Etmopterus_sheikoi")))
    phylo_list[[i]]<-temp_tree
}

write.nexus(phylo_list, file = "../../raw_data/Squaliformes_posterior_distribution.tree", translate = TRUE)

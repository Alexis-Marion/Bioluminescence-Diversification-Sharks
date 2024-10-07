library("TreeTools")
library("phytools")
library("stringr")
library("dispRity")

source("function_grafting.r")

trees<-read.tree("../../raw_data/Squaliformes_posterior_distribution.tree")

phy<-read.tree("../../raw_data/Squaliformes_extant.tree")

tab<-read.table("../../raw_data/Ts_Te_species.tsv", header = TRUE, sep ="\t")

table_trait<-read.table("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", header = TRUE, sep ="\t")

multi_tree<-function(tree, fossil_data, taxonomy, output){
phylo_list<-list()
class(phylo_list)<-"multiPhylo"
       for(i in 1:1000){
           
            temp_tree<-tree
           
            taxo_eoetmopterus_supracretaceus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Etmopterus_sheikoi", "Aculeola_nigra")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Eoetmopterus_supracretaceus", taxo_eoetmopterus_supracretaceus, 1)

            taxo_isistius_trituratus<-c("Isistius_brasiliensis")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Isistius_trituratus", taxo_isistius_trituratus, 1)

            taxo_protocentrophorus_balticus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Deania_calcea", "Centrophorus_squamosus", "Centrophorus_granulosus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protocentrophorus_balticus", taxo_protocentrophorus_balticus, 1)

            taxo_protosqualus_albertsi<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protosqualus_albertsi", taxo_protosqualus_albertsi, 1)

            taxo_protosqualus_sigei<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protosqualus_sigei", taxo_protosqualus_sigei, 1)

            taxo_squalus_minor<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Squalus_minor", taxo_squalus_minor, 1)

            taxo_trigonognathus_virginiae<-c("Trigonognathus_kabeyai")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Trigonognathus_virginiae", taxo_trigonognathus_virginiae, 1)

            taxo_eodalatias_crenulatus<-c("Dalatias_licha")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Eodalatias_crenulatus", taxo_eodalatias_crenulatus, 1)

            taxo_squalidalatias_savoiei<-extract.clade(temp_tree, getMRCA(temp_tree, c("Dalatias_licha", "Isistius_brasiliensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Squaliodalatias_savoiei", taxo_squalidalatias_savoiei, 1)

            taxo_proetmopterus_hemmooriensis<-extract.clade(temp_tree, getMRCA(temp_tree, c("Etmopterus_sheikoi", "Etmopterus_princeps", "Etmopterus_polli")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Proetmopterus_hemmooriensis", taxo_proetmopterus_hemmooriensis, 1)

            taxo_cretascymnus_westfalicus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Somniosus_microcephalus", "Rhinoscymnus_rostratus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Cretascymnus_westfalicus", taxo_cretascymnus_westfalicus, 1)

            taxo_cretascymnus_adonis<-extract.clade(temp_tree, getMRCA(temp_tree, c("Somniosus_microcephalus", "Rhinoscymnus_rostratus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Cretascymnus_adonis", taxo_cretascymnus_adonis, 1)
            phylo_list[[i]]<-temp_tree
       }
    write.nexus(phylo_list, file = output)
}

multi_tree(phy, tab, table_trait, "../../raw_data/Multi_Squaliformes_fossil_consensus_distribution.tree")

multi_tree_posterior<-function(tree, fossil_data, taxonomy, output){
phylo_list<-list()
class(phylo_list)<-"multiPhylo"
       for(i in 1:100){
           
            temp_tree<-extract.clade(tree[[i]], getMRCA(tree[[i]], c("Squalus_acanthias", "Etmopterus_princeps")))

            taxo_eoetmopterus_supracretaceus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Etmopterus_sheikoi", "Aculeola_nigra")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Eoetmopterus_supracretaceus", taxo_eoetmopterus_supracretaceus, 1)

            taxo_isistius_trituratus<-c("Isistius_brasiliensis")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Isistius_trituratus", taxo_isistius_trituratus, 1)
           
            taxo_protocentrophorus_balticus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Deania_calcea", "Centrophorus_squamosus", "Centrophorus_granulosus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protocentrophorus_balticus", taxo_protocentrophorus_balticus, 1)
            
            taxo_protosqualus_albertsi<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protosqualus_albertsi", taxo_protosqualus_albertsi, 1)

            taxo_protosqualus_sigei<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Protosqualus_sigei", taxo_protosqualus_sigei, 1)

            taxo_squalus_minor<-extract.clade(temp_tree, getMRCA(temp_tree, c("Squalus_acanthias", "Cirrhigaleus_barbifer", "Squalus_cubensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Squalus_minor", taxo_squalus_minor, 1)

            taxo_trigonognathus_virginiae<-c("Trigonognathus_kabeyai")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Trigonognathus_virginiae", taxo_trigonognathus_virginiae, 1)

            taxo_eodalatias_crenulatus<-c("Dalatias_licha")
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Eodalatias_crenulatus", taxo_eodalatias_crenulatus, 1)

            taxo_squalidalatias_savoiei<-extract.clade(temp_tree, getMRCA(temp_tree, c("Dalatias_licha", "Isistius_brasiliensis")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Squaliodalatias_savoiei", taxo_squalidalatias_savoiei, 1)

            taxo_proetmopterus_hemmooriensis<-extract.clade(temp_tree, getMRCA(temp_tree, c("Etmopterus_sheikoi", "Etmopterus_princeps", "Etmopterus_polli")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Proetmopterus_hemmooriensis", taxo_proetmopterus_hemmooriensis, 1)

            taxo_cretascymnus_westfalicus<-extract.clade(temp_tree, getMRCA(temp_tree, c("Somniosus_microcephalus", "Rhinoscymnus_rostratus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Cretascymnus_westfalicus", taxo_cretascymnus_westfalicus, 1)

            taxo_cretascymnus_adonis<-extract.clade(temp_tree, getMRCA(temp_tree, c("Somniosus_microcephalus", "Rhinoscymnus_rostratus")))$tip.label
            temp_tree<-tree_grafting_manual(temp_tree, tab, "Cretascymnus_adonis", taxo_cretascymnus_adonis, 1)
            phylo_list[[i]]<-temp_tree
       }
    write.tree(phylo_list, file = output)
}

multi_tree_posterior(trees, tab, table_trait, "../../raw_data/Multi_Squaliformes_fossil_posterior_distribution.tree")

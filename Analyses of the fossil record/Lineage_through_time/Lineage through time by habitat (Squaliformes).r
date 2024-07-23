library("ggplot2")
library("dplyr")
library("deeptime")
library("ape")
library("treeio")
library("stringr")
library("cowplot")
library("ggstream")

source("utility_function_ltt_plot.r")

phy<-read.tree("../../raw_data/Squaliformes_posterior_distribution.tree")

tab_habitat<-read.table("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep ="\t", header = TRUE)

list_files<-list.files("../Results/Genus/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE)

list_files<-list_files[grepl("se_est",list_files)]

tab_tax<-read.csv("../../raw_data/Data_taxo_genus.tsv", sep ="\t")

tab_genus_habitat<-read.table("../../raw_data/Habitat_fossil_Genus.tsv", sep ="\t", header = TRUE)

whole_tab<-tab_genus_habitat[tab_genus_habitat$Genus %in% tab_tax$Genus, ]
it<-0
for (file in list_files){
    it<- it + 1
    tab_ts_te<-read.table(file, sep ="\t", header = TRUE)
    tab_ts_te_reduced<-tab_ts_te[,-c(1,2)]
    colnames(tab_ts_te_reduced)<-c(paste(c("ts_rep_", it), collapse = ""), paste(c("te_rep_", it), collapse = ""))
    whole_tab<-cbind(whole_tab, tab_ts_te_reduced)
}

whole_tab<-whole_tab[!is.na(whole_tab$Habitat),]

write.table(whole_tab, "../../raw_data/whole_tab_genus_habitat_ts_te.tsv", sep ="\t")

tab_S <- whole_tab[whole_tab$Habitat == "S",]
tab_D <- whole_tab[whole_tab$Habitat == "D",]
tab_M <- whole_tab[whole_tab$Habitat == "M",]
tab_all <- whole_tab

trim_ts_D<-tab_D[,grepl("ts_rep_", colnames(tab_D))]
trim_te_D<-tab_D[,grepl("te_rep_", colnames(tab_D))]

trim_ts_S<-tab_S[,grepl("ts_rep_", colnames(tab_S))]
trim_te_S<-tab_S[,grepl("te_rep_", colnames(tab_S))]

trim_ts_M<-tab_M[,grepl("ts_rep_", colnames(tab_M))]
trim_te_M<-tab_M[,grepl("te_rep_", colnames(tab_M))]

trim_ts_all<-tab_all[,grepl("ts_rep_", colnames(tab_all))]
trim_te_all<-tab_all[,grepl("te_rep_", colnames(tab_all))]

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_D)){
        div<-0
        for (j in 1:nrow(trim_ts_D)){
            if (is.in.interval(trim_te_D[j,k], trim_ts_D[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_D<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_D<-as.data.frame(tab_D)
rownames(tab_D)<-NULL
colnames(tab_D)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_D_G<-tab_D

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_S)){
        div<-0
        for (j in 1:nrow(trim_ts_S)){
            if (is.in.interval(trim_te_S[j,k], trim_ts_S[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_S<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_S<-as.data.frame(tab_S)
rownames(tab_S)<-NULL
colnames(tab_S)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_S_G<-tab_S

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_M)){
        div<-0
        for (j in 1:nrow(trim_ts_M)){
            if (is.in.interval(trim_te_M[j,k], trim_ts_M[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_M<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_M<-as.data.frame(tab_M)

rownames(tab_M)<-NULL
colnames(tab_M)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_M_G<-tab_M

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_all)){
        div<-0
        for (j in 1:nrow(trim_ts_all)){
            if (is.in.interval(trim_te_all[j,k], trim_ts_all[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_all<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_all<-as.data.frame(tab_all)
rownames(tab_all)<-NULL
colnames(tab_all)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_all_G<-tab_all

list_files <- list.files("../Results/Species/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE)

list_files <- list_files[grepl("se_est",list_files)]

tab_tax_eco <- tab_habitat[,c(1:4,9)]

tab_tax_eco$Species<-gsub(" ", "_", tab_tax_eco$Species)

tab_tax<-read.csv("../../raw_data/Data_taxo_sp.tsv", sep ="\t")

tab_tax_eco <- tab_tax_eco[tab_tax_eco$Species %in% tab_tax$Species,]

whole_tab<-tab_tax_eco
it<-0
for (file in list_files){
    it<- it + 1
    tab_ts_te<-read.table(file, sep ="\t", header = TRUE)
    tab_ts_te_reduced<-tab_ts_te[,-c(1,2)]
    colnames(tab_ts_te_reduced)<-c(paste(c("ts_rep_", it), collapse = ""), paste(c("te_rep_", it), collapse = ""))
    whole_tab<-cbind(whole_tab, tab_ts_te_reduced)
}

whole_tab<-whole_tab[!is.na(whole_tab$Habitat),]

write.table(whole_tab, "../../raw_data/whole_tab_species_habitat_ts_te.tsv", sep ="\t")

tab_S<-whole_tab[whole_tab$Habitat == "S",]
tab_D<-whole_tab[whole_tab$Habitat == "D",]
tab_M<-whole_tab[whole_tab$Habitat == "M",]
tab_all <- whole_tab

trim_ts_D<-tab_D[,grepl("ts_rep_", colnames(tab_D))]
trim_te_D<-tab_D[,grepl("te_rep_", colnames(tab_D))]

trim_ts_S<-tab_S[,grepl("ts_rep_", colnames(tab_S))]
trim_te_S<-tab_S[,grepl("te_rep_", colnames(tab_S))]

trim_ts_M<-tab_M[,grepl("ts_rep_", colnames(tab_M))]
trim_te_M<-tab_M[,grepl("te_rep_", colnames(tab_M))]

trim_ts_all<-tab_all[,grepl("ts_rep_", colnames(tab_all))]
trim_te_all<-tab_all[,grepl("te_rep_", colnames(tab_all))]


tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_D)){
        div<-0
        for (j in 1:nrow(trim_ts_D)){
            if (is.in.interval(trim_te_D[j,k], trim_ts_D[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_D<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_D<-as.data.frame(tab_D)
rownames(tab_D)<-NULL
colnames(tab_D)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_D<-tab_D

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_S)){
        div<-0
        for (j in 1:nrow(trim_ts_S)){
            if (is.in.interval(trim_te_S[j,k], trim_ts_S[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_S<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_S<-as.data.frame(tab_S)
rownames(tab_S)<-NULL
colnames(tab_S)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_S<-tab_S

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_M)){
        div<-0
        for (j in 1:nrow(trim_ts_M)){
            if (is.in.interval(trim_te_M[j,k], trim_ts_M[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_M<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_M<-as.data.frame(tab_M)
rownames(tab_M)<-NULL
colnames(tab_M)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_M<-tab_M

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_all)){
        div<-0
        for (j in 1:nrow(trim_ts_all)){
            if (is.in.interval(trim_te_all[j,k], trim_ts_all[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_all<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_all<-as.data.frame(tab_all)
rownames(tab_all)<-NULL
colnames(tab_all)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_all<-tab_all

tab_habitat<-read.table("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep ="\t", header =TRUE)

list_files<-list.files("../Results/Species_combined/TS_TE_output/output_ltt/", full.names = TRUE)

list_files<-list_files[grepl("se_est",list_files)]

list_files_tax<-list.files("../Results/Species_combined/TS_TE_output/Taxonomic_equivalence/", full.names = TRUE)

colnames(tab_ts_te_reduced)<-c(paste(c("ts_rep_", it), collapse = ""), paste(c("te_rep_", it), collapse = ""))

whole_tab<-c()
for(i in 1:20){
    temp_tab_tax <- read.table(list_files_tax[i], header = TRUE, sep = "\t")
    temp_tab_ts_te <- read.table(list_files[i], header = TRUE, sep = "\t")
    temp_tab <- cbind(temp_tab_tax, temp_tab_ts_te)
    temp_tab <- temp_tab[, -c(1,3,4)]
    colnames(temp_tab)<-c("Species", paste(c("ts_rep_", i), collapse = ""), paste(c("te_rep_", i), collapse = ""))
    if(length(whole_tab) == 0){
        whole_tab <- temp_tab
    }
    else{
        whole_tab <- merge(whole_tab, temp_tab, by = "Species")
    }
}

ecological_cat<-c()
for( i in 1:nrow(tab_habitat)){
    temp_sp_name <- gsub(" ", "_", tab_habitat$Species[i])
    if( temp_sp_name %in% whole_tab$Species){
        temp_vec<-c(tab_habitat$Order[i], tab_habitat$Family[i], tab_habitat$Genus[i], temp_sp_name, tab_habitat$Habitat[i])
        ecological_cat <- rbind(ecological_cat, temp_vec)
    }
}
ecological_cat <- as.data.frame(ecological_cat)
colnames(ecological_cat) <- c("Order", "Family", "Genus", "Species", "Habitat")

data_combined <- merge(ecological_cat, whole_tab, by = "Species")
data_combined <- cbind(data_combined[c(2,3,4)], data_combined[c(1)], data_combined[-c(1,2,3,4)])
data_combined <- data_combined[!is.na(data_combined$Habitat),]

tab_S<-data_combined[data_combined$Habitat == "S",]
tab_D<-data_combined[data_combined$Habitat == "D",]
tab_M<-data_combined[data_combined$Habitat == "M",]
tab_all<-data_combined

trim_ts_D<-tab_D[,grepl("ts_rep_", colnames(tab_D))]
trim_te_D<-tab_D[,grepl("te_rep_", colnames(tab_D))]
trim_ts_D<-apply(trim_ts_D, 2, as.numeric)
trim_te_D<-apply(trim_te_D, 2, as.numeric)

trim_ts_S<-tab_S[,grepl("ts_rep_", colnames(tab_S))]
trim_te_S<-tab_S[,grepl("te_rep_", colnames(tab_S))]
trim_ts_S<-apply(trim_ts_S, 2, as.numeric)
trim_te_S<-apply(trim_te_S, 2, as.numeric)

trim_ts_M<-tab_M[,grepl("ts_rep_", colnames(tab_M))]
trim_te_M<-tab_M[,grepl("te_rep_", colnames(tab_M))]
trim_ts_M<-apply(trim_ts_M, 2, as.numeric)
trim_te_M<-apply(trim_te_M, 2, as.numeric)

trim_ts_all<-tab_all[,grepl("ts_rep_", colnames(tab_all))]
trim_te_all<-tab_all[,grepl("te_rep_", colnames(tab_all))]
trim_ts_all<-apply(trim_ts_all, 2, as.numeric)
trim_te_all<-apply(trim_te_all, 2, as.numeric)

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_D)){
        div<-0
        for (j in 1:nrow(trim_ts_D)){
            if (is.in.interval(trim_te_D[j,k], trim_ts_D[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_D<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_D<-as.data.frame(tab_D)
rownames(tab_D)<-NULL
colnames(tab_D)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_combined_D<-tab_D

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_S)){
        div<-0
        for (j in 1:nrow(trim_ts_S)){
            if (is.in.interval(trim_te_S[j,k], trim_ts_S[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_S<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_S<-as.data.frame(tab_S)
rownames(tab_S)<-NULL
colnames(tab_S)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_combined_S<-tab_S

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_M)){
        div<-0
        for (j in 1:nrow(trim_ts_M)){
            if (is.in.interval(trim_te_M[j,k], trim_ts_M[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_M<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_M<-as.data.frame(tab_M)
rownames(tab_M)<-NULL
colnames(tab_M)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_combined_M<-tab_M

tab<-c()
for(i in seq(-145, 0, by = 0.5)){
    div_vec<-c(i)
    for (k in 1:ncol(trim_ts_all)){
        div<-0
        for (j in 1:nrow(trim_ts_all)){
            if (is.in.interval(trim_te_all[j,k], trim_ts_all[j,k], i, i + 0.5)){
                div <- div + 1      
            }
        }
        div_vec<-c(div_vec, div)
    }
    tab<-rbind(tab, div_vec)
}

tab_all<-cbind(-rev(tab[,1]),
          rev(apply(tab, 1, median)),
          rev(apply(tab, 1, quantile_ts_te)[2,]),
          rev(apply(tab, 1, quantile_ts_te)[1,]))

tab_all<-as.data.frame(tab_all)
rownames(tab_all)<-NULL
colnames(tab_all)<-c("Age","Diversity", "min_Diversity", "max_Diversity")

ltt_combined_all<-tab_all

write.table(ltt_all_G, "ltt_gn.tsv", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ltt_all, "ltt_sp.tsv", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(ltt_combined_all, "ltt_sp_combined.tsv", sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

tab <- cbind((ltt_D_G/ltt_all_G)[,2], (ltt_D/ltt_all)[,2], (ltt_combined_D/ltt_combined_all)[,2])

write.table(tab, "data_prop.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

data_1<-data_maker(ltt_D, ltt_S, ltt_M)
data_1$Time_bin<- (-data_1$Time_bin)
data_2<-data_maker(ltt_combined_D, ltt_combined_S, ltt_combined_M)
data_2$Time_bin<- (-data_2$Time_bin)
data_3<-data_maker(ltt_D_G, ltt_S_G, ltt_M_G)
data_3$Time_bin<- (-data_3$Time_bin)

b1<-ggplot(data_1, aes(fill=Habitat, y=Sample, x=Time_bin)) + 
    geom_bar(position="fill", stat="identity", width = 5.1) + 
scale_fill_manual("legend", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species",
       x = "Year",
       y = "Proportion", 
       fill = "Habitat")

b2<-ggplot(data_2, aes(fill=Habitat, y=Sample, x=Time_bin)) + 
    geom_bar(position="fill", stat="identity", width = 5.1) + 
scale_fill_manual("legend", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species combined",
       x = "Year",
       y = "Proportion", 
       fill = "Habitat")

b3<-ggplot(data_3, aes(fill=Habitat, y=Sample, x=Time_bin)) + 
    geom_bar(position="fill", stat="identity", width = 5.1) + 
scale_fill_manual("legend", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) +
  labs(title = "Habitat occupation through time",
       subtitle = "Genus",
       x = "Year",
       y = "Proportion", 
       fill = "Habitat")

p1<-ggplot(data_1, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="ridge", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

p3<-ggplot(data_2, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="ridge", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species combined",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

p5<-ggplot(data_3, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="ridge", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Genus",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

p2<-ggplot(data_1, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="proportional", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

p4<-ggplot(data_2, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="proportional", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Species combined",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

p6<-ggplot(data_3, aes(x = Time_bin,
                     y = Sample, 
                     fill = Habitat)) +
  geom_stream(type="proportional", n_grid = 100000, true_range = c("both"), extra_span = 0) +
  labs(title = "Habitat occupation through time",
       subtitle = "Genus",
       x = "Year",
       y = "Diversity (nb. lineages)", 
       fill = "Habitat") + 
scale_fill_manual("Habitat", values = c("D" = "#5E4fA2", "S" = "#74BDCA", "M" = "#3A8B9B")) + 
  theme_minimal() 

pdf("all_habitat_ltt.pdf", width = 30, height = 10)
    plot_grid(b1, p1, p2, b2, p3, p4, b3, p5, p6, ncol = 3, rel_widths = c(1/3, 1/3, 1/3))
dev.off()

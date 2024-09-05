library("ggplot2")
library("dplyr")
library("deeptime")
library("ape")
library("treeio")
library("stringr")
library("cowplot")
library("ggstream")

source("utility_function_ltt_plot.r")

data_whole_gen<-seq(0, 145, 0.5)
data_whole_sp<-seq(0, 145, 0.5)
data_whole_sp_comb<-seq(0, 145, 0.5)

for(l in 1:1000){
### Loading genus data

list_files<-list.files("../Results/Genus/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE)
    
tab_habitat<-read.table("../../raw_data/Trait_data_Squaliformes_Fossil.tsv", sep ="\t", header =TRUE)

list_files<-list_files[grepl("se_est",list_files)]

tab_tax<-read.csv("../../raw_data/Data_taxo_genus.tsv", sep ="\t")

tab_genus_habitat<-read.table("../raw_data/Habitat_fossil_Genus.tsv", sep ="\t", header = TRUE)

### Preparing dataframe with shared Ts/Te

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

### Lineage lifespan by habitat
whole_tab$Habitat<-sample(whole_tab$Habitat, length(whole_tab$Habitat))
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

## Species

### Loading species data

list_files <- list.files("../Results/Species/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE)

list_files <- list_files[grepl("se_est",list_files)]

tab_tax_eco <- tab_habitat[,c(1:4,9)]

tab_tax_eco$Species<-gsub(" ", "_", tab_tax_eco$Species)

tab_tax<-read.csv("../../raw_data/Data_taxo_species.tsv", sep ="\t")

tab_tax_eco <- tab_tax_eco[tab_tax_eco$Species %in% tab_tax$Species,]

### Preparing dataframe with shared Ts/Te

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

### Lineage lifespan by habitat (Fossil only)

whole_tab$Habitat<-sample(whole_tab$Habitat, length(whole_tab$Habitat))
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

### Lineage lifespan by habitat (Combined)

data_combined <- merge(ecological_cat, whole_tab, by = "Species")
data_combined <- cbind(data_combined[c(2,3,4)], data_combined[c(1)], data_combined[-c(1,2,3,4)])
data_combined <- data_combined[!is.na(data_combined$Habitat),]

data_combined$Habitat<-sample(data_combined$Habitat, length(data_combined$Habitat))
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
    
data_whole_gen<-cbind(data_whole_gen, propmaker(ltt_D_G$Diversity, ltt_S_G$Diversity, ltt_M_G$Diversity))
data_whole_sp<-cbind(data_whole_sp, propmaker(ltt_D$Diversity, ltt_S$Diversity, ltt_M$Diversity))
data_whole_sp_comb<-cbind(data_whole_sp_comb, propmaker(ltt_combined_D$Diversity, ltt_combined_S$Diversity, ltt_combined_M$Diversity))
}

write.table(data_whole_gen, "data_whole_gen.tsv", quote = FALSE, sep ="\t", row.names =FALSE)
write.table(data_whole_sp, "data_whole_sp.tsv", quote = FALSE, sep ="\t", row.names =FALSE)
write.table(data_whole_sp_comb, "data_whole_sp_comb.tsv", quote = FALSE,  sep ="\t", row.names =FALSE)

tab_prop<-read.table("data_prop.tsv", sep ="\t", header =TRUE)

tab<-apply(data_whole_gen[,-1], 1, quantile, probs = c(0.025,0.975))
vec<-seq(0, 145, by = 0.5)
vec_T_F<-rep("No Difference", 291)
for(i in 1:291){
    if(round(tab_prop[i,1],3) > round(tab[2,i],3)){
        vec_T_F[i]<-("Upper") 
    }
    if(round(tab_prop[i,1],3) < round(tab[1,i],3)){
        vec_T_F[i]<-("Lower") 
    }
}
df_gen<-as.data.frame(cbind(vec, round(t(tab),3),round(tab_prop[,1],3), vec_T_F))
colnames(df_gen)<-c("Age", "2.5%", "97.5%", "Empirical proportion", "Significance")

tab<-apply(data_whole_sp[,-1], 1, quantile, probs = c(0.025,0.975))
vec<-seq(0, 145, by = 0.5)
vec_T_F<-rep("No Difference", 291)
for(i in 1:291){
    if(round(tab_prop[i,2],3) > round(tab[2,i],3)){
        vec_T_F[i]<-("Upper") 
    }
    if(round(tab_prop[i,2],3) < round(tab[1,i],3)){
        vec_T_F[i]<-("Lower") 
    }
}
df_sp<-as.data.frame(cbind(vec, round(t(tab),3),round(tab_prop[,2],3), vec_T_F))
colnames(df_sp)<-c("Age", "2.5%", "97.5%", "Empirical proportion", "Significance")

tab<-apply(data_whole_sp_comb[,-1], 1, quantile, probs = c(0.025,0.975))
vec<-seq(0, 145, by = 0.5)
vec_T_F<-rep("No Difference", 291)
for(i in 1:291){
    if(round(tab_prop[i,3],3) > round(tab[2,i],3)){
        vec_T_F[i]<-("Upper") 
    }
    if(round(tab_prop[i,3],3) < round(tab[1,i],3)){
        vec_T_F[i]<-("Lower") 
    }
}
df_sp_comb<-as.data.frame(cbind(vec, round(t(tab),3),round(tab_prop[,3],3), vec_T_F))
colnames(df_sp_comb)<-c("Age", "2.5%", "97.5%", "Empirical proportion", "Significance")

write.table(df_gen, "df_gen.tsv", row.names = FALSE, quote = FALSE, sep ="\t")
write.table(df_sp, "df_sp.tsv", row.names = FALSE, quote = FALSE, sep ="\t")
write.table(df_sp_comb, "df_sp_comb.tsv", row.names = FALSE, quote = FALSE, sep ="\t")

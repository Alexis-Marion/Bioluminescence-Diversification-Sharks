library("stringr")
library("phytools")
library("dispRity")
library("dplyr")
library("deeptime")
library("ggplot2")
library("treeio")
library("stringr")
library("cowplot")

phy_list <- read.tree("../../raw_data/Squaliformes_posterior_distribution.tree")

table_habitat_sp <- read.table("../../raw_data/whole_tab_species_ts_te.tsv", sep ="\t")

table_habitat_gn <- read.table("../../raw_data/whole_tab_genus_ts_te.tsv", sep ="\t")

Squalidae <-c('Squalus_blainvillei', 'Squalus_megalops', 'Squalus_albicaudus', 'Squalus_raoulensis', 'Squalus_formosus', 'Squalus_crassispinus', 'Squalus_brevirostris', 'Squalus_hemipinnis', 'Squalus_edmundsi', 'Squalus_mitsukurii', 'Squalus_japonicus', 'Squalus_nasutus', 'Squalus_grahami', 'Squalus_cubensis', 'Squalus_montalbani', 'Squalus_chloroculus', 'Squalus_acanthias', 'Squalus_suckleyi', 'Cirrhigaleus_australis', 'Cirrhigaleus_barbifer', 'Cirrhigaleus_asper')
Centrophoridae <-c('Centrophorus_squamosus', 'Centrophorus_granulosus', 'Centrophorus_zeehaani', 'Centrophorus_uyato', 'Centrophorus_atromarginatus', 'Centrophorus_longipinnis', 'Centrophorus_isodon', 'Centrophorus_tessellatus', 'Centrophorus_westraliensis', 'Centrophorus_harrissoni', 'Centrophorus_lesliei', 'Centrophorus_lusitanicus', 'Deania_hystricosa', 'Deania_calcea', 'Deania_profundorum', 'Deania_quadrispinosa')
Dalatiidae <-c('Squaliolus_aliae', 'Euprotomicrus_bispinatus', 'Squaliolus_laticaudus', 'Dalatias_licha', 'Isistius_brasiliensis')
Somniosidae <-c('Centroscymnus_owstonii', 'Oxynotus_paradoxus', 'Oxynotus_bruniensis', 'Oxynotus_centrina', 'Scymnodon_ringens', 'Centroscymnus_coelolepis', 'Scymnodon_ichiharai', 'Scymnodalatias_albicauda', 'Centroselachus_crepidater', 'Zameus_squamulosus', 'Somniosus_microcephalus', 'Rhinoscymnus_rostratus')
Etmopteridae <-c('Etmopterus_brachyurus', 'Etmopterus_samadiae', 'Etmopterus_sheikoi', 'Etmopterus_evansi', 'Etmopterus_brosei', 'Etmopterus_molleri', 'Etmopterus_pycnolepis', 'Etmopterus_lucifer', 'Etmopterus_dislineatus', 'Etmopterus_alphus', 'Etmopterus_bullisi', 'Etmopterus_bigelowi', 'Etmopterus_pusillus', 'Etmopterus_fusus', 'Etmopterus_splendidus', 'Etmopterus_pseudosqualiolus', 'Etmopterus_sentosus', 'Etmopterus_unicolor', 'Etmopterus_princeps', 'Etmopterus_spinax', 'Etmopterus_viator', 'Etmopterus_litvinovi', 'Etmopterus_granulosus', 'Etmopterus_compagnoi', 'Etmopterus_dianthus', 'Etmopterus_schultzi', 'Etmopterus_polli', 'Etmopterus_virens', 'Etmopterus_gracilispinis', 'Etmopterus_carteri', 'Centroscyllium_fabricii', 'Centroscyllium_ritteri', 'Centroscyllium_granulatum', 'Centroscyllium_nigrum', 'Aculeola_nigra', 'Trigonognathus_kabeyai')
Etmopterus <-c('Etmopterus_brachyurus','Etmopterus_samadiae','Etmopterus_sheikoi','Etmopterus_evansi','Etmopterus_brosei','Etmopterus_molleri','Etmopterus_pycnolepis','Etmopterus_lucifer','Etmopterus_dislineatus','Etmopterus_alphus','Etmopterus_bullisi','Etmopterus_bigelowi','Etmopterus_pusillus','Etmopterus_fusus','Etmopterus_splendidus','Etmopterus_pseudosqualiolus','Etmopterus_sentosus','Etmopterus_unicolor','Etmopterus_princeps','Etmopterus_spinax','Etmopterus_viator','Etmopterus_litvinovi','Etmopterus_granulosus','Etmopterus_compagnoi','Etmopterus_dianthus','Etmopterus_schultzi','Etmopterus_polli','Etmopterus_virens','Etmopterus_gracilispinis','Etmopterus_carteri')
Somniosidae_reduced <-c('Centroscymnus_owstonii','Oxynotus_paradoxus','Oxynotus_bruniensis','Oxynotus_centrina','Scymnodon_ringens','Centroscymnus_coelolepis','Scymnodon_ichiharai','Scymnodalatias_albicauda','Centroselachus_crepidater','Zameus_squamulosus')
Bioluminescent <- c(Dalatiidae, Somniosidae, Etmopteridae)

table_habitat_maker <- function(tab){
    for (i in 1:nrow(tab)){
        temp_tab <- tab
        vec_habitat <- c()
        vec_age <- c()
        if("D" %in% tab$Habitat){
            temp_tab_D <- temp_tab %>% filter(Habitat == "D")
            temp_tab_D <- temp_tab_D[,grepl("ts_rep_", colnames(temp_tab_D))]
            temp_tab_D <- temp_tab_D[which.max(apply(temp_tab_D, 1, FUN=mean)),]
            vec_habitat <- c("D")
            vec_age <- c(t(temp_tab_D))
        }
        if("M" %in% tab$Habitat){
            temp_tab_M <- temp_tab %>% filter(Habitat == "M")
            temp_tab_M <- temp_tab_M[,grepl("ts_rep_", colnames(temp_tab_M))]
            temp_tab_M <- temp_tab_M[which.max(apply(temp_tab_M, 1, FUN=mean)),]
            vec_habitat <- c(vec_habitat, "M")
            vec_age <- c(vec_age, t(temp_tab_M))
        }
        if("S" %in% tab$Habitat){
            temp_tab_S <- temp_tab %>% filter(Habitat == "S")
            temp_tab_S <- temp_tab_S[,grepl("ts_rep_", colnames(temp_tab_S))]
            temp_tab_S <- temp_tab_S[which.max(apply(temp_tab_S, 1, FUN=mean)),]
            vec_habitat <- c(vec_habitat, "S")
            vec_age <- c(vec_age, t(temp_tab_S))
        }
        return(data.frame(habitat=factor(rep(vec_habitat, each=20)), age = vec_age))
    }
}

table_biol_maker <- function (tab, phy_list, taxa_list){
    TO_C <- c()
    TO_S <- c()
    for (i in 1:100) {
        phy <- phy_list[[i]]
        tree_ages <- tree.age(phy, order = "past", fossil = TRUE, 
            digits = 3)
        if(length(taxa_list) == 1){
            TO_C <- c(TO_C, 0)
            TO_S <- c(TO_S, tree_ages[tree_ages[, 2] == phy$edge[phy$edge[,2] == which(phy$tip.label == taxa_list), 1], ][, 1])
        }
        else{
            TO_C <- c(TO_C, tree_ages[tree_ages[, 2] == getMRCA(phy, taxa_list), ][, 1])
            TO_S <- c(TO_S, tree_ages[tree_ages[, 2] == phy$edge[phy$edge[,2] == getMRCA(phy, taxa_list), 1], ][, 1])
            
        }
    }
    temp_tab <- tab
    temp_tab <- temp_tab[, grepl("ts_rep_", colnames(temp_tab))]
    temp_tab <- temp_tab[which.max(apply(temp_tab, 1, FUN = mean)), ]
    
    return(data.frame(habitat = factor(rep(c("TO_S", "TO_S", "TO_S", "TO_S", "TO_S", "F", "TO_C", "TO_C", "TO_C", "TO_C", "TO_C"), 
        each = 20)), age = c(TO_S, t(temp_tab),  TO_C)))
}

table_deep_biol_maker <- function (tab, phy_list, taxa_list){
    TO_C <- c()
    TO_S <- c()
    for (i in 1:100) {
        phy <- phy_list[[i]]
        tree_ages <- tree.age(phy, order = "past", fossil = TRUE, 
            digits = 3)
        if(length(taxa_list) == 1){
            TO_C <- c(TO_C, 0)
            TO_S <- c(TO_S, tree_ages[tree_ages[, 2] == phy$edge[phy$edge[,2] == which(phy$tip.label == taxa_list), 1], ][, 1])
        }
        else{
            TO_C <- c(TO_C, tree_ages[tree_ages[, 2] == getMRCA(phy, taxa_list), ][, 1])
            TO_S <- c(TO_S, tree_ages[tree_ages[, 2] == phy$edge[phy$edge[,2] == getMRCA(phy, taxa_list), 1], ][, 1])
            
        }
    }
    temp_tab <- tab
    temp_tab_D <- temp_tab %>% filter(Habitat == "D")
    temp_tab_D <- temp_tab_D[,grepl("ts_rep_", colnames(temp_tab_D))]
    temp_tab_D <- temp_tab_D[which.max(apply(temp_tab_D, 1, FUN=mean)),]
    
    return(data.frame(habitat = factor(rep(c("TO_S", "TO_S", "TO_S", "TO_S", "TO_S", "F", "TO_C", "TO_C", "TO_C", "TO_C", "TO_C"), 
        each = 20)), age = c(TO_S, t(temp_tab_D),  TO_C)))
}

plot_maker_age <- function(tab, min_age, max_age){
    p<-ggplot(tab, aes(x=age, fill=habitat)) +
        geom_density(alpha=0.4) + 
        xlim(max_age, min_age) +
        theme_bw() +
        theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) + 
        scale_fill_manual(values=c("#D2B48C", "#405BBC", "#40826D")) + 
        coord_geo(dat = list("stage", "epochs"), pos = as.list(rep("bottom", 2)), abbrv = FALSE, height = list(unit(1, "lines"), unit(1, "lines")),  size = list(2.5, 2.5))
    return(p)
}

plot_maker_age_biol <- function(tab, min_age, max_age){
    p<-ggplot(tab, aes(x=age, fill=habitat)) +
        geom_density(alpha=0.4) + 
        xlim(max_age, min_age) +
        theme_bw() +
        theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))  + 
        scale_fill_manual(values=c("#D2B48C", "#40826D", "#405BBC")) +
        coord_geo(dat = list("stage", "epochs"), pos = as.list(rep("bottom", 2)), abbrv = FALSE, height = list(unit(1, "lines"), unit(1, "lines")),  size = list(2.5, 2.5))
        return(p)
}

table_habitat_sp_Etmopteridae <- table_habitat_sp %>% filter(Family == "Etmopteridae")

table_habitat_sp_Somniosidae <- table_habitat_sp %>% filter(Family == "Somniosidae")

table_habitat_sp_Dalatiidae <- table_habitat_sp %>% filter(Family == "Dalatiidae")

table_habitat_sp_Centrophoridae <- table_habitat_sp %>% filter(Family == "Centrophoridae")

table_habitat_sp_Squalidae <- table_habitat_sp %>% filter(Family == "Squalidae")

table_habitat_sp_Zameus <- table_habitat_sp %>% filter(Genus == "Zameus")

table_habitat_sp_Bioluminescent <- table_habitat_sp %>% filter(Family != "Centrophoridae") %>% filter(Family != "Squalidae")

tab_Squalidae_sp <- table_habitat_maker(table_habitat_sp_Squalidae)

tab_Centrophoridae_sp <- table_habitat_maker(table_habitat_sp_Centrophoridae)

tab_Dalatiidae_sp <- table_habitat_maker(table_habitat_sp_Dalatiidae)

tab_Somniosidae_sp <- table_habitat_maker(table_habitat_sp_Somniosidae)

tab_Etmopteridae_sp <- table_habitat_maker(table_habitat_sp_Etmopteridae)

biol_Squalidae_sp <- table_biol_maker(table_habitat_sp_Squalidae, phy_list, Squalidae)

biol_Centrophoridae_sp <- table_biol_maker(table_habitat_sp_Centrophoridae, phy_list, Centrophoridae)

biol_Dalatiidae_sp <- table_biol_maker(table_habitat_sp_Dalatiidae, phy_list, Dalatiidae)

biol_Somniosidae_sp <- table_biol_maker(table_habitat_sp_Somniosidae, phy_list, Somniosidae)

biol_Etmopteridae_sp <- table_biol_maker(table_habitat_sp_Etmopteridae, phy_list, Etmopteridae)

biol_Zameus_sp <- table_biol_maker(table_habitat_sp_Zameus, phy_list, "Zameus_squamulosus")

deep_biol_Bioluminescent_sp <- table_deep_biol_maker(table_habitat_sp_Bioluminescent, phy_list, Bioluminescent)

deep_biol_Dalatiidae_sp <- table_deep_biol_maker(table_habitat_sp_Dalatiidae, phy_list, Dalatiidae)

deep_biol_Etmopteridae_sp <- table_deep_biol_maker(table_habitat_sp_Etmopteridae, phy_list, Etmopteridae)

deep_biol_Zameus_sp <- table_deep_biol_maker(table_habitat_sp_Zameus, phy_list, "Zameus_squamulosus")

b1_s<- plot_maker_age(tab_Squalidae_sp, 66, 145)

b2_s<- plot_maker_age(tab_Centrophoridae_sp, 56, 100.5)

b3_s<- plot_maker_age(tab_Dalatiidae_sp, 56, 100.5)

b4_s<- plot_maker_age(tab_Somniosidae_sp, 66, 100.5)

b5_s<- plot_maker_age(tab_Etmopteridae_sp, 66, 100.5)

p1_s<-plot_maker_age_biol(biol_Squalidae_sp, 20, 160)

p2_s<-plot_maker_age_biol(biol_Centrophoridae_sp, 30, 160)

p3_s<-plot_maker_age_biol(biol_Dalatiidae_sp, 66, 160)

p4_s<-plot_maker_age_biol(biol_Somniosidae_sp, 80, 150)

p5_s<-plot_maker_age_biol(biol_Etmopteridae_sp, 66, 150)

p6_s<-plot_maker_age_biol(biol_Zameus_sp, 0, 66)

bp1_s<-plot_maker_age_biol(deep_biol_Bioluminescent_sp, 75, 180)

bp2_s<-plot_maker_age_biol(deep_biol_Dalatiidae_sp, 66, 160)

bp3_s<-plot_maker_age_biol(deep_biol_Etmopteridae_sp, 66, 160)

bp4_s<-plot_maker_age_biol(deep_biol_Zameus_sp, 0, 66)

table_habitat_gn_Etmopteridae <- table_habitat_gn %>% filter(Family == "Etmopteridae")

table_habitat_gn_Somniosidae <- table_habitat_gn %>% filter(Family == "Somniosidae")

table_habitat_gn_Dalatiidae <- table_habitat_gn %>% filter(Family == "Dalatiidae")

table_habitat_gn_Centrophoridae <- table_habitat_gn %>% filter(Family == "Centrophoridae")

table_habitat_gn_Squalidae <- table_habitat_gn %>% filter(Family == "Squalidae")

table_habitat_gn_Zameus <- table_habitat_gn %>% filter(Genus == "Zameus")

table_habitat_gn_Bioluminescent <- table_habitat_gn %>% filter(Family != "Centrophoridae") %>% filter(Family != "Squalidae")

tab_Squalidae_gn <- table_habitat_maker(table_habitat_gn_Squalidae)

tab_Centrophoridae_gn <- table_habitat_maker(table_habitat_gn_Centrophoridae)

tab_Dalatiidae_gn <- table_habitat_maker(table_habitat_gn_Dalatiidae)

tab_Somniosidae_gn <- table_habitat_maker(table_habitat_gn_Somniosidae)

tab_Etmopteridae_gn <- table_habitat_maker(table_habitat_gn_Etmopteridae)

biol_Squalidae_gn <- table_biol_maker(table_habitat_gn_Squalidae, phy_list, Squalidae)

biol_Centrophoridae_gn <- table_biol_maker(table_habitat_gn_Centrophoridae, phy_list, Centrophoridae)

biol_Dalatiidae_gn <- table_biol_maker(table_habitat_gn_Dalatiidae, phy_list, Dalatiidae)

biol_Somniosidae_gn <- table_biol_maker(table_habitat_gn_Somniosidae, phy_list, Somniosidae)

biol_Etmopteridae_gn <- table_biol_maker(table_habitat_gn_Etmopteridae, phy_list, Etmopteridae)

biol_Zameus_gn <- table_biol_maker(table_habitat_gn_Zameus, phy_list, "Zameus_squamulosus")

deep_biol_Bioluminescent_gn <- table_deep_biol_maker(table_habitat_gn_Bioluminescent, phy_list, Bioluminescent)

deep_biol_Dalatiidae_gn <- table_deep_biol_maker(table_habitat_gn_Dalatiidae, phy_list, Dalatiidae)

deep_biol_Etmopteridae_gn <- table_deep_biol_maker(table_habitat_gn_Etmopteridae, phy_list, Etmopteridae)

deep_biol_Zameus_gn <- table_deep_biol_maker(table_habitat_gn_Zameus, phy_list, "Zameus_squamulosus")

b1_g <- plot_maker_age(tab_Squalidae_gn, 66, 145)

b2_g <- plot_maker_age(tab_Centrophoridae_gn, 56, 100.5)

b3_g <- plot_maker_age(tab_Dalatiidae_gn, 56, 100.5)

b4_g <- plot_maker_age(tab_Somniosidae_gn, 66, 100.5)

b5_g <- plot_maker_age(tab_Etmopteridae_gn, 66, 100.5)

p1_g <- plot_maker_age_biol(biol_Squalidae_gn, 20, 160)

p2_g <- plot_maker_age_biol(biol_Centrophoridae_gn, 30, 160)

p3_g <- plot_maker_age_biol(biol_Dalatiidae_gn, 66, 160)

p4_g <- plot_maker_age_biol(biol_Somniosidae_gn, 80, 150)

p5_g <- plot_maker_age_biol(biol_Etmopteridae_gn, 66, 150)

p6_g <- plot_maker_age_biol(biol_Zameus_gn, 0, 66)

bp1_g<-plot_maker_age_biol(deep_biol_Bioluminescent_gn, 75, 160)

bp2_g<-plot_maker_age_biol(deep_biol_Dalatiidae_gn, 66, 160)

bp3_g<-plot_maker_age_biol(deep_biol_Etmopteridae_gn, 66, 160)

bp4_g<-plot_maker_age_biol(deep_biol_Zameus_gn, 0, 66)

pdf("tempo_deep_shallow.pdf", width = 10, height = 20)
    plot_grid(b1_s, b1_g, b2_s, b2_g, b3_s, b3_g, b4_s, b4_g, b5_s, b5_g, ncol = 2, rel_widths = c(1/2, 1/2))
dev.off()

pdf("tempo_origination.pdf", width = 10, height = 20)
    plot_grid(p1_s, p1_g, p2_s, p2_g, p3_s, p3_g, p4_s, p4_g, p5_s, p5_g, p6_s, p6_g, ncol = 2, rel_widths = c(1/2, 1/2))
dev.off()

pdf("tempo_biol_deep.pdf", width = 10, height = 20)
    plot_grid(bp1_s, bp1_g, bp2_s, bp2_g, bp3_s, bp3_g, bp4_s, bp4_g, ncol = 2, rel_widths = c(1/2, 1/2))
dev.off()

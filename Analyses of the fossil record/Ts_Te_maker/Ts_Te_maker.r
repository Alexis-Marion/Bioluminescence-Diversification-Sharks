data_taxo_species <- read.table("../../raw_data/Data_taxo_species.tsv", header = TRUE, sep ="\t")

data_taxo_genus <- read.table("../../raw_data/Data_taxo_genus.tsv", header = TRUE, sep ="\t")

list_files_species <- list.files("../Results/Species/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE, pattern = "se_est")

list_files_genus <- list.files("../Results/Genus/5_MA_bins_EXT_EP/pyrate_mcmc_logs/", full.names = TRUE, pattern = "se_est")

vec_Ts_Te <- c()
for(i in 1:20){
    vec_Ts_Te<- c(vec_Ts_Te, (paste("Ts_rep_",i, sep = "")), (paste("Te_rep_",i, sep = "")))
}

Ts <- c()
Te <- c()
Ts_Te <- data_taxo_species

for(file in list_files_species){
    file <- read.table(file, header = TRUE)
    Ts <- cbind(Ts, file[,3])
    Te <- cbind(Te, file[,4])
    Ts_Te <- cbind(Ts_Te, file[,3:4])
}

tab_sp <- cbind(data_taxo_species,
apply(X = Ts, 1, FUN = mean),
apply(X = Te, 1, FUN = mean),
apply(X = Ts, 1, FUN = sd),
apply(X = Te, 1, FUN = sd))
colnames(tab_sp)[5:8] <- c("Ts", "Te", "Ts_sd", "Te_sd")

whole_tab_sp <- Ts_Te
colnames(whole_tab_sp)[5:44] <- vec_Ts_Te

write.table(tab_sp, "../../raw_data/TS_TE_SP_replicated.tsv", sep ="\t", quote = FALSE, row.names = FALSE)

write.table(whole_tab_sp, "../../raw_data/whole_tab_species_ts_te.tsv", sep ="\t", quote = FALSE, row.names = FALSE)

Ts_Te <- data_taxo_genus

for(file in list_files_genus){
    file <- read.table(file, header = TRUE)
    Ts_Te <- cbind(Ts_Te, file[,3:4])
}

whole_tab_gn <- Ts_Te
colnames(whole_tab_gn)[4:43] <- vec_Ts_Te

write.table(whole_tab_sp, "../../raw_data/whole_tab_genus_ts_te.tsv", sep ="\t", quote = FALSE, row.names = FALSE)

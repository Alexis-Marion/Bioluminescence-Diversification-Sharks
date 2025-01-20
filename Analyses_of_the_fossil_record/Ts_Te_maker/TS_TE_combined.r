library("ape")

tree_list<-read.tree("../../raw_data/Squaliformes_posterior_distribution.tree")

tabl<-read.table("Species_combined/Occ_sp_Squali_TaxonList.txt", sep ="\t", header = TRUE)

lst_files<-list.files("Species_combined/TS_TE/", full.names = TRUE)

tr_nb_list<-sample(100, 20, replace = FALSE)

for (i in 1:20){
    tr<-tree_list[[tr_nb_list[i]]]
    temp<-read.table(lst_files[[i]], sep ="\t", header = TRUE)
    temp_tabl<-cbind(tabl, temp)
    temp_phy<-as.data.frame(cbind(tr$tip.label,tr$edge.length[tr$edge[,2] <= length(tr$tip.label)]))
    temp_phy2<- temp_phy[! temp_phy[,1] %in% temp_tabl[,1],2]
    temp_neonto<- as.data.frame(cbind(rep(0,length(temp_phy2)), ((nrow(temp_tabl)+1):(nrow(temp_tabl)+length(temp_phy2))), temp_phy2, rep(0,length(temp_phy2))))
    colnames(temp_neonto)<-colnames(temp)
    write.table(rbind(temp, temp_neonto), paste("Species_combined/TS_TE_output/Occ_sp_Squali_combined", i, "_G_KEEP_BDS_se_est.txt", sep =""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

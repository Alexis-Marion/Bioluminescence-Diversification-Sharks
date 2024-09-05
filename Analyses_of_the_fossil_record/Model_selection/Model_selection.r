library("dplyr")

args = commandArgs(trailingOnly=TRUE)

list_files_1<-list.files(args[1], full.names = TRUE)

list_files_1<-list_files_1[grepl("KEEP_BDS_mcmc.log",list_files_1)]

list_files_2<-list.files(args[2], full.names = TRUE)

list_files_2<-list_files_2[grepl("KEEP_BDS_mcmc.log",list_files_2)]

compute<-function(list_file_1, list_file_2, k, len, funct){
data_df<-c()
    for (i in 1:len){
        
        var_conv_1 = "FALSE"
        var_conv_2 = "FALSE"
        
        file_1<-(list_file_1[grepl(paste("_", i, "_", sep = ""), list_file_1)])
        file_2<-(list_file_2[grepl(paste("_", i, "_", sep = ""), list_file_2)])
        
        if(grepl("KEEP", file_1)){
            var_conv_1 <- "TRUE"
        }
        
        if(grepl("KEEP", file_1)){
            var_conv_2 <- "TRUE"
        }
        
        file_1<-read.table(file_1, sep ="\t", header = TRUE)
        file_2<-read.table(file_2, sep ="\t", header = TRUE)
        
        file_1<-file_1[-c(1:nrow(file_1)*k),]

        l1_tot<-file_1[,grepl("lambda_",colnames(file_1))]
        m1_tot<-file_1[,grepl("mu_",colnames(file_1))]
        TS1_tot<-file_1[,grepl("_TS",colnames(file_1))]
        TE1_tot<-file_1[,grepl("_TE",colnames(file_1))]
        
        file_2<-file_2[-c(1:nrow(file_2)*k),]

        l2_tot<-file_2[,grepl("lambda_",colnames(file_2))]
        m2_tot<-file_2[,grepl("mu_",colnames(file_2))]
        TS2_tot<-file_2[,grepl("_TS",colnames(file_2))]
        TE2_tot<-file_2[,grepl("_TE",colnames(file_2))]
        
        if(ncol(TS1_tot) != ncol(TS2_tot)){
            
            maxlist_TS<-list(TS1_tot,TS2_tot)
            maxlist_TE<-list(TE1_tot,TE2_tot)
            
            TSmax<-maxlist_TS[[which.max(c(ncol(TS1_tot),ncol(TS2_tot)))]]
            TSmin<-maxlist_TS[[which.min(c(ncol(TS1_tot),ncol(TS2_tot)))]]
            TEmax<-maxlist_TE[[which.max(c(ncol(TE1_tot),ncol(TE2_tot)))]]
            TEmin<-maxlist_TE[[which.min(c(ncol(TE1_tot),ncol(TE2_tot)))]]
            
            maxlist_TS[[which.max(c(ncol(TS1_tot),ncol(TS2_tot)))]]<-TSmax[,colnames(TSmax) %in% colnames(TSmin)]
            maxlist_TE[[which.max(c(ncol(TE1_tot),ncol(TE2_tot)))]]<-TEmax[,colnames(TEmax) %in% colnames(TEmin)]
            
            TS1_tot<-maxlist_TS[[1]]
            TS2_tot<-maxlist_TS[[2]]
            TE1_tot<-maxlist_TE[[1]]
            TE2_tot<-maxlist_TE[[2]]
        }
        
        lst<-list(list(l1_tot, l2_tot), list(m1_tot, m2_tot), list(TS1_tot, TS2_tot), list(TE1_tot, TE2_tot))
        
        for (j in 1:4){
            if (j == 1){
                var <- "Lambda"
            }
            if (j == 2){
                var <- "Mu"
            }
            if (j == 3){
                var <- "TS"
            }
            if (j == 4){
                var <- "TE"
            }
            temp_data1<-apply(X = lst[[j]][[1]], MARGIN = 2, FUN = funct)
            temp_data2<-apply(X = lst[[j]][[2]], MARGIN = 2, FUN = funct)
            if(shapiro.test(temp_data1)$p<0.05|shapiro.test(temp_data2)$p<0.05){
                mean_diff<-wilcox.test
                test<-"Wilcoxon"
            }else{
                mean_diff<-t.test
                test<-"T.test"
            }
            res <- mean_diff(temp_data1, temp_data2, paired = TRUE, conf.int = TRUE, conf.level = 0.95, exact = TRUE)
            data_df<-rbind(data_df,c(i, var, test, round(res$statistic, 3), round(res$p.value, 3), round(res$estimate, 3), round(res$conf.int[1], 3), round(res$conf.int[2], 3), var_conv_1, var_conv_2))
        }
    }
    data_df<-as.data.frame(data_df)
    colnames(data_df)<-c("Replicates", "Estimated Parameters", "Test", "Statistic", "p.value","Mean/Median difference", "Lower interval", "Upper interval", "Conv_1", "Conv_2")
    return(data_df)
}

write.table(compute(list_files_1, list_files_2, 0.01, 20, median), args[3], sep ="\t", row.names = FALSE)

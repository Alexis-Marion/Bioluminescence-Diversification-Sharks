mround <- function(x,base){
    base*round(x/base)
}
quantile_ts_te<-function(data_ts_te){
        return(c(sort(data_ts_te)[0.975*length(data_ts_te)], sort(data_ts_te, decreasing = TRUE)[0.975*length(data_ts_te)]))
}
is.in.interval<-function(age1, age2, interval1, interval2){
    temp_vec<-(-seq(round(age1, 1), round(age2, 1), by = 0.1))
    if(interval1 %in% temp_vec | interval2 %in% temp_vec){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}
get_age_bound2<-function(data){
mn<-as.numeric(rownames(deeptime::epochs[deeptime::epochs[,3]<=min(data$Age),][nrow(deeptime::epochs[deeptime::epochs[,3]<=min(data$Age),]),]))
mx<-as.numeric(rownames(deeptime::epochs[deeptime::epochs[,3]<=max(data$Age),][nrow(deeptime::epochs[deeptime::epochs[,3]<=max(data$Age),]),]))
  
output<-list(deeptime::epochs[mn:mx,2],deeptime::epochs[mn:mx,])
    return(output)
}

ltt_plot <- function(ltt_df, #has to be in the format returned by the `extract_ltt()` function from the `extract_param_from_PyRate_outputs.R` script
                     x_breaks = get_age_bound2(ltt_df)[[1]],
                     y_breaks = seq(0,20,5),
                     y_labels = seq(0,20,5),
                     y_limits = c(0, 20),
                     main=NA,
                     x_lab = "Time (Ma)",
                     y_lab = "Diversity (nb. lineages)",
                     geoscale = get_age_bound2(ltt_df)[[2]],
                     abbr = TRUE){
  temp_data<-data.frame(cbind(get_age_bound2(ltt_df)[[1]]), rep(1, length (get_age_bound2(ltt_df)[[1]])))
  colnames(temp_data)<-c("X", "Y")
  # Proper plot
  p <- ggplot(data = ltt_df, aes(x = Age, y = Diversity)) +
    scale_x_reverse(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks,
                       limits = y_limits) +
    geom_line(data = temp_data, aes(x = X, y = Y), linewidth = 1, colour = "#FFFFFF", linetype ="dashed") +
    geom_ribbon(aes(x = Age, ymin = min_Diversity, ymax = max_Diversity), 
                fill = "#74BDCA",
                alpha = 0.8) +
    geom_line(linewidth = 1, colour = "#3A8B9B") +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    coord_geo(dat = geoscale, abbrv = abbr, size = 5)
    return(p)
}

ltt_plot_two_states <- function(ltt_df, ltt_df1, ltt_df2, #has to be in the format returned by the `extract_ltt()` function from the `extract_param_from_PyRate_outputs.R` script
                     x_breaks = get_age_bound2(ltt_df)[[1]],
                     y_breaks = seq(0,15,5),
                     y_labels = seq(0,15,5),
                     y_limits = c(0, 15),
                     main=NA,
                     x_lab = "Time (Ma)",
                     y_lab = "Diversity (nb. lineages)",
                     geoscale = get_age_bound2(ltt_df)[[2]],
                     abbr = TRUE){
  temp_data<-data.frame(cbind(get_age_bound2(ltt_df)[[1]]), rep(0.5, length (get_age_bound2(ltt_df)[[1]])))
  colnames(temp_data)<-c("X", "Y")
  # Proper plot
  p <- ggplot(data = ltt_df, aes(x = Age, y = Diversity)) +
    scale_x_reverse(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks,
                       limits = y_limits) +
    geom_line(data = temp_data, aes(x = X, y = Y), linewidth = 1, colour = "#FFFFFF", linetype ="dashed") +
    geom_ribbon(data = ltt_df1, aes(x = Age, ymin = min_Diversity, ymax = max_Diversity), 
                fill = "#CBC3E3",
                alpha = 0.8) +
    geom_line(data = ltt_df1, aes(x = Age), linewidth = 1, colour = "#5E4fA2") +
        geom_ribbon(data = ltt_df2, aes(x = Age, ymin = min_Diversity, ymax = max_Diversity), 
                fill = "#CCEAFD",
                alpha = 0.8) +
    geom_line(data = ltt_df2, aes(x = Age), linewidth = 1, colour = "#74BDCB") + 
    xlab(x_lab) +
    ylab(y_lab) +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    coord_geo(dat = geoscale, abbrv = abbr, size = 5)
    return(p)
}

ltt_plot_ratio <- function(ltt_df, ltt_df_ratio, #has to be in the format returned by the `extract_ltt()` function from the `extract_param_from_PyRate_outputs.R` script
                     x_breaks = get_age_bound2(ltt_df)[[1]],
                     y_breaks = seq(0,10,1),
                     y_labels = seq(0,10,1),
                     y_limits = c(0, 10),
                     main=NA,
                     x_lab = "Time (Ma)",
                     y_lab = "Ratio of shallow-to-deep diversity",
                     geoscale = get_age_bound2(ltt_df)[[2]],
                     abbr = TRUE){
  temp_data<-data.frame(cbind(get_age_bound2(ltt_df)[[1]]), rep(1, length (get_age_bound2(ltt_df)[[1]])))
  colnames(temp_data)<-c("X", "Y")
  # Proper plot
  p <- ggplot(data = ltt_df, aes(x = Age, y = Diversity)) +
    scale_x_reverse(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks,
                       limits = y_limits) +
    geom_line(data = temp_data, aes(x = X, y = Y), linewidth = 1, colour = "#A9A9A9", linetype ="dashed") +
    geom_ribbon(data = ltt_df_ratio, aes(x = Age, ymin = min_Diversity, ymax = max_Diversity), 
                fill = "#C4B8C8",
                alpha = 0.8) +
    geom_line(data = ltt_df_ratio, aes(x = Age),linewidth = 1, colour = "#979AAA") +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    coord_geo(dat = geoscale, abbrv = abbr, size = 5)
    return(p)
}

combined_ltt<-function(ltt_extant, ltt_extinct){
    if (which.max(c(nrow(ltt_extant), nrow(ltt_extinct))) == 1){
        ltt_total <- ltt_extant
        ltt_minus <- ltt_extinct
    }
    else{
        ltt_total <- ltt_extinct
        ltt_minus <- ltt_extant
    }
    for (i in 1:nrow(ltt_minus)){
        ltt_total[i,c(2,3,4)] <- ltt_minus[i,c(2,3,4)] + ltt_total[i,c(2,3,4)]
    }
    return(ltt_total)
}

data_maker<-function(ltt_D, ltt_S, ltt_M){
data<-c()
for(i in 1:nrow(ltt_D)){
    vec_D<-c(ltt_D[i,1], "D", ltt_D[i,2])
    if(ltt_D[i,1] %in% ltt_S[i,1]){
        vec_S<-c(ltt_D[i,1], "S", ltt_S[i,2])  
    }
    else{
        vec_S<-c(ltt_D[i,1], "S", 0) 
    }
    if(ltt_D[i,1] %in% ltt_M[i,1]){
        vec_M<-c(ltt_D[i,1], "M", ltt_M[i,2])
    }
    else{
        vec_M<-c(ltt_D[i,1], "M", 0) 
    }
data<-rbind(data,vec_D, vec_S, vec_M)
data<-as.data.frame(data)
colnames(data)<-c("Time_bin", "Habitat", "Sample")
data$Sample<-as.numeric(data$Sample)
data$Time_bin<-(as.numeric(data$Time_bin))   
} 
return(data)
} 

propmaker<-function(vecA, vecB, vecC){
    vecTot<- (vecA + vecB + vecC)
    prop<-vecTot
    prop[!vecTot==0]<- vecA[!vecTot==0]/vecTot[!vecTot==0]
    return(prop)
}

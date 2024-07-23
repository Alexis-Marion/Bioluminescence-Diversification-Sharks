

pdf(file='../Results/Genus/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_gn_Squali_6_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_6_G_KEEP_BDS_mcmc"
ts=c(46.60924125654578, 49.15729895720128,65.48503225978008,86.97384513660195,70.17607721507389,86.53841887757667,38.46491625787677,88.58989833603393,86.03662334825924,58.21527351510295,66.78461631993721,36.05829109562769,67.0842662128209,85.66058014246029,48.96697305490782,47.406550433536346,47.13340418953835,23.258832506859232,76.1166445429621,86.28491792135864,59.33395504172703,67.60939275290131,69.46533551953206,55.68564461463383,33.63506409788083,49.70724308814847,77.68104329131533,97.00874719582359,130.20775341852453,94.10654255459148,83.75920397143062,49.7311181731226,36.232182376796125,42.76103771614775,96.73791374932735,49.01788827304256,99.86681005008796,49.284293972233705,38.45286784768543)
te=c(40.49587655255478, 11.904070996732107,60.643078461773854,65.53876851862663,0.0,0.0,0.0,68.645014181322,65.60889109721752,0.0,0.0,16.089484457703243,43.973866837411975,67.3995604125956,13.135731366126913,0.0,0.0,0.0,70.25172996989164,61.95794349970116,0.0,49.02948731889377,67.39553226527583,11.632498607469616,0.0,12.05848716527896,64.95048369983924,65.32112266895096,67.24627017220668,66.97146569077617,0.0,0.0,0.0,0.0,13.040000874799423,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,14,15,14,15,14,13,12,13,12,13,12,11,12,11,10,9,8,9,10,11,12,13,14,15,14,15,16,17,18,19,18,19,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(130.20775341852453, 99.86681005008796,97.00874719582359,96.73791374932735,94.10654255459148,88.58989833603393,86.97384513660195,86.53841887757667,86.28491792135864,86.03662334825924,85.66058014246029,83.75920397143062,77.68104329131533,76.1166445429621,70.25172996989164,70.17607721507389,69.46533551953206,68.645014181322,67.60939275290131,67.3995604125956,67.39553226527583,67.24627017220668,67.0842662128209,66.97146569077617,66.78461631993721,65.60889109721752,65.53876851862663,65.48503225978008,65.32112266895096,64.95048369983924,61.95794349970116,60.643078461773854,59.33395504172703,58.21527351510295,55.68564461463383,49.7311181731226,49.70724308814847,49.284293972233705,49.15729895720128,49.02948731889377,49.01788827304256,48.96697305490782,47.406550433536346,47.13340418953835,46.60924125654578,43.973866837411975,42.76103771614775,40.49587655255478,38.46491625787677,38.45286784768543,36.232182376796125,36.05829109562769,33.63506409788083,23.258832506859232,16.089484457703243,13.135731366126913,13.040000874799423,12.05848716527896,11.904070996732107,11.632498607469616,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
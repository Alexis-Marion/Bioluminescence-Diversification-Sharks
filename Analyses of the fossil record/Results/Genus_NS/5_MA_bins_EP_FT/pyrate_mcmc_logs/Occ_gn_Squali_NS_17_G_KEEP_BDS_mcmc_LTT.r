

pdf(file='../Results/Genus_NS/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_NS_17_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_17_G_KEEP_BDS_mcmc"
ts=c(50.07074085597268, 67.04414029851876,89.17659821236431,74.39397902909135,86.57873757608154,35.605151988021774,88.11075631523524,85.31206404597528,59.4905656764725,69.43923697382185,31.984047122332992,67.07115116980314,85.49839957173977,50.483887967211146,47.41633599679827,87.58189496462566,59.10977360137458,70.22045043175301,57.35090867239495,35.93287811138267,50.155503336023386,76.23564371610877,96.45682625966211,130.75969628083436,94.11308989779518,83.66303642110189,50.1807684571838,37.331774867592806,43.44245129127713,96.40954979100354,47.71907089515889,98.68340129276145,50.84991203050926,32.57923398230468)
te=c(9.589605127151268, 60.103550755698045,65.97382222657633,0.0,0.0,0.0,68.2690146894105,65.53514800194993,0.0,0.0,18.98687902549451,44.228142404548635,66.09929527244584,12.424942643185952,0.0,58.55492587243646,0.0,50.449013418411845,10.451189221742718,0.0,12.137408776042603,66.54184107117699,66.7015487356209,68.67071137953359,66.34707424097863,0.0,0.0,0.0,0.0,16.13056454551949,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,15,14,15,16,15,14,13,12,11,10,9,10,11,10,11,12,13,12,13,14,15,16,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(130.75969628083436, 98.68340129276145,96.45682625966211,96.40954979100354,94.11308989779518,89.17659821236431,88.11075631523524,87.58189496462566,86.57873757608154,85.49839957173977,85.31206404597528,83.66303642110189,76.23564371610877,74.39397902909135,70.22045043175301,69.43923697382185,68.67071137953359,68.2690146894105,67.07115116980314,67.04414029851876,66.7015487356209,66.54184107117699,66.34707424097863,66.09929527244584,65.97382222657633,65.53514800194993,60.103550755698045,59.4905656764725,59.10977360137458,58.55492587243646,57.35090867239495,50.84991203050926,50.483887967211146,50.449013418411845,50.1807684571838,50.155503336023386,50.07074085597268,47.71907089515889,47.41633599679827,44.228142404548635,43.44245129127713,37.331774867592806,35.93287811138267,35.605151988021774,32.57923398230468,31.984047122332992,18.98687902549451,16.13056454551949,12.424942643185952,12.137408776042603,10.451189221742718,9.589605127151268,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
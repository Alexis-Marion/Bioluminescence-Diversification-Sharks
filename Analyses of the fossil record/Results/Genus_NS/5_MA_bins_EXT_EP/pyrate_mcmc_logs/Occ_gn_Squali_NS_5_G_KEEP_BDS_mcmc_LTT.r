

pdf(file='../Results/Genus_NS/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_gn_Squali_NS_5_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_5_G_KEEP_BDS_mcmc"
ts=c(49.5388916075922, 66.81779412920079,88.74678559408854,67.297342978051,85.57721697833617,39.05743172316753,88.63976535701897,85.5480692538789,60.581377441856404,63.98242348298687,35.358091391875924,64.3264906439335,85.80181627726031,49.698867362265766,47.334139255448335,88.44137665038356,58.65311886907234,70.85496951090923,56.22464416559627,35.964455995397046,49.35570794768311,74.46863011449175,96.63114235970164,129.48464875020326,94.23805995928959,85.32912305258202,50.424615503338536,35.9536140984683,43.394288527882686,97.11798649862062,47.92188810654794,96.89610711609814,50.08677481431116,35.51934546128384)
te=c(12.503793456352316, 60.924222925135425,67.49648483989965,0.0,0.0,0.0,69.1494000296359,65.16372893881388,0.0,0.0,16.91323714712848,44.303495310939816,65.89085343936877,15.50320580432261,0.0,61.7905319967104,0.0,46.9118779010542,12.816116751688432,0.0,11.821103709376468,65.31734384642431,69.58156624183326,68.91539427079903,67.20487754267907,0.0,0.0,0.0,0.0,13.792861137406938,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,12,11,10,11,10,11,10,9,8,9,10,9,8,9,10,11,12,13,14,15,16,17,18,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(129.48464875020326, 97.11798649862062,96.89610711609814,96.63114235970164,94.23805995928959,88.74678559408854,88.63976535701897,88.44137665038356,85.80181627726031,85.57721697833617,85.5480692538789,85.32912305258202,74.46863011449175,70.85496951090923,69.58156624183326,69.1494000296359,68.91539427079903,67.49648483989965,67.297342978051,67.20487754267907,66.81779412920079,65.89085343936877,65.31734384642431,65.16372893881388,64.3264906439335,63.98242348298687,61.7905319967104,60.924222925135425,60.581377441856404,58.65311886907234,56.22464416559627,50.424615503338536,50.08677481431116,49.698867362265766,49.5388916075922,49.35570794768311,47.92188810654794,47.334139255448335,46.9118779010542,44.303495310939816,43.394288527882686,39.05743172316753,35.964455995397046,35.9536140984683,35.51934546128384,35.358091391875924,16.91323714712848,15.50320580432261,13.792861137406938,12.816116751688432,12.503793456352316,11.821103709376468,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
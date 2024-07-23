

pdf(file='../Results/Genus_NS/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_gn_Squali_NS_8_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_8_G_KEEP_BDS_mcmc"
ts=c(49.83067795101515, 65.53669042964593,90.71983232682125,73.61654366765025,87.31279595313481,38.89832790713927,88.87051504071306,85.80406905216218,62.386385077327766,63.619940989350766,35.944981736279736,64.02825437547735,85.24443198154168,49.622141430465255,47.20263138372371,88.32129513417195,58.60720028558856,68.39486509886618,60.44086051498734,36.33301060590962,49.4283639951131,75.12671804082001,97.7682350733119,128.17529182973016,94.23793433424615,83.83729679596743,49.99776488416056,34.40938390857173,41.43619248853706,97.60637656070301,48.2672682208247,100.6892453907907,49.58618694425612,33.776446308836235)
te=c(10.866583832644467, 61.44090361697263,65.59944325245743,0.0,0.0,0.0,68.98845738266202,65.3149193989989,0.0,0.0,16.990048269503234,43.84472397257913,67.92155222420095,12.939465495601087,0.0,61.70828084644153,0.0,50.43092281486673,10.610650595712563,0.0,11.686134148579953,65.35048892772821,67.5659639310836,69.1978968222073,66.4085135734244,0.0,0.0,0.0,0.0,12.892154643243856,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,12,13,12,11,10,9,10,9,8,9,10,11,10,9,10,11,10,11,12,13,14,15,16,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(128.17529182973016, 100.6892453907907,97.7682350733119,97.60637656070301,94.23793433424615,90.71983232682125,88.87051504071306,88.32129513417195,87.31279595313481,85.80406905216218,85.24443198154168,83.83729679596743,75.12671804082001,73.61654366765025,69.1978968222073,68.98845738266202,68.39486509886618,67.92155222420095,67.5659639310836,66.4085135734244,65.59944325245743,65.53669042964593,65.35048892772821,65.3149193989989,64.02825437547735,63.619940989350766,62.386385077327766,61.70828084644153,61.44090361697263,60.44086051498734,58.60720028558856,50.43092281486673,49.99776488416056,49.83067795101515,49.622141430465255,49.58618694425612,49.4283639951131,48.2672682208247,47.20263138372371,43.84472397257913,41.43619248853706,38.89832790713927,36.33301060590962,35.944981736279736,34.40938390857173,33.776446308836235,16.990048269503234,12.939465495601087,12.892154643243856,11.686134148579953,10.866583832644467,10.610650595712563,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
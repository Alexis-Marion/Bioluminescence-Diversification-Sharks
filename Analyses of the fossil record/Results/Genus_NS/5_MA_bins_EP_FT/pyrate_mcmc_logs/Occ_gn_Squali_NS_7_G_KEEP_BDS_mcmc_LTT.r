

pdf(file='../Results/Genus_NS/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_NS_7_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_7_G_KEEP_BDS_mcmc"
ts=c(50.19308480378447, 67.53296241506953,89.75879000075825,68.99943913215473,86.78043723754531,36.16916121039272,88.66025616265418,86.2527141004214,60.01278739589402,65.83597750813236,32.97074060647539,66.49520089903042,85.87356564104512,49.754058877788125,47.51337855226394,87.03175141748507,58.70432654868442,67.55655527547425,56.39007148909878,36.44710774858029,50.224965307624174,78.27841223613476,96.21593230423814,129.12357255842673,93.21979299491215,84.12621467716795,49.40240991943568,35.436967582828906,44.989454049325545,96.48697174825787,49.65905056708356,97.72290094233672,49.9657139161534,33.66371608943541)
te=c(12.247434083197728, 60.86726267205147,67.05854807621947,0.0,0.0,0.0,68.75294237320914,65.68032710512387,0.0,0.0,17.265228829997916,46.9480699356305,64.8737310820014,14.842166907174457,0.0,57.75517924113167,0.0,50.095868527342574,11.937736917882953,0.0,11.67616429344279,66.20478565029374,68.2016234010968,69.17680783993566,66.85894970984275,0.0,0.0,0.0,0.0,13.991123136112808,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,12,13,12,11,12,13,12,11,12,11,12,11,10,9,10,11,10,11,12,13,12,13,14,15,16,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(129.12357255842673, 97.72290094233672,96.48697174825787,96.21593230423814,93.21979299491215,89.75879000075825,88.66025616265418,87.03175141748507,86.78043723754531,86.2527141004214,85.87356564104512,84.12621467716795,78.27841223613476,69.17680783993566,68.99943913215473,68.75294237320914,68.2016234010968,67.55655527547425,67.53296241506953,67.05854807621947,66.85894970984275,66.49520089903042,66.20478565029374,65.83597750813236,65.68032710512387,64.8737310820014,60.86726267205147,60.01278739589402,58.70432654868442,57.75517924113167,56.39007148909878,50.224965307624174,50.19308480378447,50.095868527342574,49.9657139161534,49.754058877788125,49.65905056708356,49.40240991943568,47.51337855226394,46.9480699356305,44.989454049325545,36.44710774858029,36.16916121039272,35.436967582828906,33.66371608943541,32.97074060647539,17.265228829997916,14.842166907174457,13.991123136112808,12.247434083197728,11.937736917882953,11.67616429344279,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                


pdf(file='../Results/Genus_NS/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_NS_14_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_14_G_KEEP_BDS_mcmc"
ts=c(48.55917791826654, 66.6812916068787,87.38153411909092,68.49152035557528,86.13537324120253,37.292781138979436,89.06230545904906,85.44674813719827,62.403090119627734,66.16457295879111,32.42939888488741,65.45628896332785,85.25843439541643,48.00267702263974,47.23403635482511,87.0819037927029,59.894358480648116,68.31895242185028,52.54394432424686,35.59251289048253,48.94980572349299,78.08245316813976,96.86029842151257,130.9261922668265,94.05787027296999,84.22458931539924,49.37237204470038,37.05449604499729,43.96617577242415,96.80755505392067,47.249387288691494,98.20066627185611,49.18127035682218,31.14659996597537)
te=c(12.647693670386955, 62.02200906296197,65.53457194842315,0.0,0.0,0.0,67.95700569327762,65.25744656257415,0.0,0.0,19.23287572982304,45.94500672263363,65.71932900337946,13.399811145523506,0.0,58.02571880977827,0.0,48.43476592002401,12.640568230436088,0.0,12.035610252023428,65.24550771697255,68.08216669128304,69.01621702142356,65.72586697856148,0.0,0.0,0.0,0.0,15.342978162700526,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,12,13,14,13,12,13,14,13,12,11,12,11,10,11,10,11,10,11,12,13,14,15,14,15,16,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(130.9261922668265, 98.20066627185611,96.86029842151257,96.80755505392067,94.05787027296999,89.06230545904906,87.38153411909092,87.0819037927029,86.13537324120253,85.44674813719827,85.25843439541643,84.22458931539924,78.08245316813976,69.01621702142356,68.49152035557528,68.31895242185028,68.08216669128304,67.95700569327762,66.6812916068787,66.16457295879111,65.72586697856148,65.71932900337946,65.53457194842315,65.45628896332785,65.25744656257415,65.24550771697255,62.403090119627734,62.02200906296197,59.894358480648116,58.02571880977827,52.54394432424686,49.37237204470038,49.18127035682218,48.94980572349299,48.55917791826654,48.43476592002401,48.00267702263974,47.249387288691494,47.23403635482511,45.94500672263363,43.96617577242415,37.292781138979436,37.05449604499729,35.59251289048253,32.42939888488741,31.14659996597537,19.23287572982304,15.342978162700526,13.399811145523506,12.647693670386955,12.640568230436088,12.035610252023428,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                


pdf(file='../Results/Genus_NS/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_gn_Squali_NS_14_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_14_G_KEEP_BDS_mcmc"
ts=c(48.42471296308391, 66.79168366301084,87.38556946348291,68.35504143626859,86.03490113086062,38.08017820499389,88.78342488945509,85.3828476311542,62.290581116790044,66.2158055707562,33.57542216594975,65.24649843692804,85.37259174733373,49.29236902550305,47.10191687594819,86.90114719582435,59.6990426605096,68.3442492459556,52.27447736750333,36.39289729165242,49.04128164488137,78.1852999015133,96.79185207423717,131.01916167877653,94.15950988335021,84.68032685004601,49.416600241078385,37.19685561099922,44.561314111445895,96.76380632842542,47.14797445555383,98.22120057734861,49.39972160716209,30.964981224395995)
te=c(12.650046555664531, 62.11274560796694,65.38772639281609,0.0,0.0,0.0,68.41570021233807,65.17378989284995,0.0,0.0,19.100362234647413,45.47520611043132,65.36317158801997,13.620110084256567,0.0,59.41027780045194,0.0,48.45022986444838,12.77584675908838,0.0,12.07537775641157,65.0686420212866,68.58766045162221,69.27159066940939,65.30412731070179,0.0,0.0,0.0,0.0,15.26713809023359,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,12,11,10,11,12,13,14,13,12,11,12,11,10,11,10,11,10,11,12,13,14,15,14,15,16,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(131.01916167877653, 98.22120057734861,96.79185207423717,96.76380632842542,94.15950988335021,88.78342488945509,87.38556946348291,86.90114719582435,86.03490113086062,85.3828476311542,85.37259174733373,84.68032685004601,78.1852999015133,69.27159066940939,68.58766045162221,68.41570021233807,68.35504143626859,68.3442492459556,66.79168366301084,66.2158055707562,65.38772639281609,65.36317158801997,65.30412731070179,65.24649843692804,65.17378989284995,65.0686420212866,62.290581116790044,62.11274560796694,59.6990426605096,59.41027780045194,52.27447736750333,49.416600241078385,49.39972160716209,49.29236902550305,49.04128164488137,48.45022986444838,48.42471296308391,47.14797445555383,47.10191687594819,45.47520611043132,44.561314111445895,38.08017820499389,37.19685561099922,36.39289729165242,33.57542216594975,30.964981224395995,19.100362234647413,15.26713809023359,13.620110084256567,12.77584675908838,12.650046555664531,12.07537775641157,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
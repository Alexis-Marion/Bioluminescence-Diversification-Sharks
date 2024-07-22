

pdf(file='../Results/Genus/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_20_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_20_G_KEEP_BDS_mcmc"
ts=c(46.56204128707723, 49.5238337640267,66.9645773158868,90.03035267070516,68.46922202654619,87.04276747533248,40.87189209499043,89.25218141298568,86.22296331904472,59.90095739594679,64.98631538422269,35.48517759237381,66.10748656739699,86.22287788144496,51.33435188454702,47.470802341349064,47.117653567111574,25.013930496741924,78.07068446058311,87.15504104221986,59.47508067250766,68.12698542055169,74.29375873389124,58.20631717099439,36.059555916414595,50.365217388692564,78.10181737054296,97.15662864114535,131.55794224205957,93.82321514168179,80.48020084970949,50.092901118973835,35.07484429912634,40.92352876097879,96.71839024333381,50.24871533155274,98.64427329672311,50.48172637854474,34.09351303359941)
te=c(40.470964543463886, 9.535394004344111,63.07270217579374,65.91164603957886,0.0,0.0,0.0,69.48275702227275,65.53538674273778,0.0,0.0,19.189599571296124,44.80295288386158,65.760487898002,12.800608881791552,0.0,0.0,0.0,73.40419691724689,61.098405388076905,0.0,48.35912153795564,69.47968355873077,11.730864624073646,0.0,10.617780769090473,65.46854366810307,68.84333510171919,66.06073272523713,67.2506200024635,0.0,0.0,0.0,0.0,13.394439044016709,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,13,12,11,12,13,12,13,14,13,12,11,10,9,10,9,8,9,10,11,12,13,14,15,16,17,16,17,18,19,18,19,20,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(131.55794224205957, 98.64427329672311,97.15662864114535,96.71839024333381,93.82321514168179,90.03035267070516,89.25218141298568,87.15504104221986,87.04276747533248,86.22296331904472,86.22287788144496,80.48020084970949,78.10181737054296,78.07068446058311,74.29375873389124,73.40419691724689,69.48275702227275,69.47968355873077,68.84333510171919,68.46922202654619,68.12698542055169,67.2506200024635,66.9645773158868,66.10748656739699,66.06073272523713,65.91164603957886,65.760487898002,65.53538674273778,65.46854366810307,64.98631538422269,63.07270217579374,61.098405388076905,59.90095739594679,59.47508067250766,58.20631717099439,51.33435188454702,50.48172637854474,50.365217388692564,50.24871533155274,50.092901118973835,49.5238337640267,48.35912153795564,47.470802341349064,47.117653567111574,46.56204128707723,44.80295288386158,40.92352876097879,40.87189209499043,40.470964543463886,36.059555916414595,35.48517759237381,35.07484429912634,34.09351303359941,25.013930496741924,19.189599571296124,13.394439044016709,12.800608881791552,11.730864624073646,10.617780769090473,9.535394004344111,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
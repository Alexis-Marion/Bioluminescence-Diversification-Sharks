

pdf(file='../Results/Genus/5_MA_bins_EP/pyrate_mcmc_logs/Occ_gn_Squali_6_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_6_G_KEEP_BDS_mcmc"
ts=c(46.66010702474124, 49.72054326639324,65.17017602110533,87.00474919050916,70.03735088530047,86.61575310486023,37.4794628732708,88.05048579214652,85.97250742258618,58.20689401351815,66.18253877676055,35.54190974082901,67.0294350240007,85.70885462040197,48.973578894257756,47.21098380368319,47.15152618278456,23.264411153201117,75.8741603343992,86.45772433405543,59.28843670145431,67.61546268856976,69.4282587577588,56.37768882775695,39.73689686705999,49.72352402842786,77.936504477663,96.75515628545085,130.26722623210586,94.05866905172896,83.9401819509995,49.218042947320214,34.82059260023022,42.492997050096946,96.70277897825792,49.01973131447295,99.8902173226092,48.90517948068538,35.17337230967027)
te=c(40.115386709271746, 11.879614104819355,60.48702890252486,65.64434396541145,0.0,0.0,0.0,69.28525718826262,65.72544911088507,0.0,0.0,16.389291071277828,45.32976638135858,67.46245457545548,13.197547991205091,0.0,0.0,0.0,70.27571367229707,61.49315844909782,0.0,48.85804927669124,67.47387490399693,11.483465458016857,0.0,12.075819113378342,65.15456044014851,65.81786176753464,67.31170133598299,67.29281307066243,0.0,0.0,0.0,0.0,13.116092437680729,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,14,15,14,15,14,13,12,11,12,13,12,11,10,11,10,9,8,9,10,11,12,13,14,15,16,17,16,17,18,19,18,19,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(130.26722623210586, 99.8902173226092,96.75515628545085,96.70277897825792,94.05866905172896,88.05048579214652,87.00474919050916,86.61575310486023,86.45772433405543,85.97250742258618,85.70885462040197,83.9401819509995,77.936504477663,75.8741603343992,70.27571367229707,70.03735088530047,69.4282587577588,69.28525718826262,67.61546268856976,67.47387490399693,67.46245457545548,67.31170133598299,67.29281307066243,67.0294350240007,66.18253877676055,65.81786176753464,65.72544911088507,65.64434396541145,65.17017602110533,65.15456044014851,61.49315844909782,60.48702890252486,59.28843670145431,58.20689401351815,56.37768882775695,49.72352402842786,49.72054326639324,49.218042947320214,49.01973131447295,48.973578894257756,48.90517948068538,48.85804927669124,47.21098380368319,47.15152618278456,46.66010702474124,45.32976638135858,42.492997050096946,40.115386709271746,39.73689686705999,37.4794628732708,35.54190974082901,35.17337230967027,34.82059260023022,23.264411153201117,16.389291071277828,13.197547991205091,13.116092437680729,12.075819113378342,11.879614104819355,11.483465458016857,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
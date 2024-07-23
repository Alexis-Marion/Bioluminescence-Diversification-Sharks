

pdf(file='../Results/Species_NS/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_sp_Squali_NS_7_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_7_G_KEEP_BDS_mcmc"
ts=c(46.24403708467408, 66.52690337068246,65.7251022318572,84.68969064354873,86.78374424839049,21.988532919974254,29.21126778122316,7.123469518316583,83.86887029925393,78.57377420170891,7.141557023777552,70.32927827695119,83.55593967359819,19.74634033109611,20.30849345643626,28.489974551963993,53.560738865212485,59.065610442186255,69.85942307711323,84.362693482687,46.42957608331301,29.022850642484276,69.14086375332474,20.980566874180656,58.80309769494216,66.12230108383231,32.77984613475874,7.016493166431487,22.394614234240507,48.513973787720296,77.9169464450095,67.78378341300639,78.0841337171392,99.12001933899704,102.13374395578131,109.92875779843575,111.6086756972543,90.57907160633387,5.829948710825167,63.78363614146185,49.49481237650922,46.3615058994251,30.726208949497632,74.47720432824183,7.331909962305469,66.78525863074725,34.04816268562471,71.51964410518605,69.55147536261082,65.91757385914839,19.3880316352086,1.542223433042018,82.25913785642317,54.39680124169527,54.25245650872249,77.90894306307297,21.456182591648155)
te=c(37.949833321222364, 61.45296728507076,62.30449110997799,65.16378772956423,82.59452997800469,0.0,16.33066901298038,0.0,60.73368334408425,66.09131883238062,0.0,65.5572881219508,70.18793173762184,0.0,0.0,20.794801286381936,46.41066289369482,48.30298145651841,67.24810324476097,69.0141158548792,40.19733654499189,14.109011849235593,64.72368724181246,2.8331078875374702,34.49296714297766,48.58382578787084,27.29334572747111,0.0,10.878236981995371,39.34772794320706,66.77975577116197,65.6433511104564,73.89586723204289,94.73967931933974,94.4860920073734,103.93127570212403,94.80139662782356,69.99087107813214,0.0,60.48864151073824,38.08280966130497,38.16488146339577,12.597884421701972,69.6456277345401,0.0,52.68822141847846,13.096867445263703,61.01337476048234,65.20051293168005,39.08760694159392,1.5900106632198174,0.0,69.25393416090424,39.393425714136725,39.12261331212584,70.2258017750965,0.0)
div_traj=c(0, 1,2,1,2,3,2,1,0,1,2,3,4,5,6,5,6,7,8,9,10,11,10,11,12,11,10,9,10,9,10,9,10,9,10,9,10,9,10,11,10,11,12,11,10,9,8,7,8,7,6,5,4,3,4,5,6,7,8,7,8,7,8,7,8,7,8,9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,5,6,7,8,9,8,9,10,11,10,9,8,7,6,7,8,9,10,11,10,9,10)
time_events=c(111.6086756972543, 109.92875779843575,103.93127570212403,102.13374395578131,99.12001933899704,94.80139662782356,94.73967931933974,94.4860920073734,90.57907160633387,86.78374424839049,84.68969064354873,84.362693482687,83.86887029925393,83.55593967359819,82.59452997800469,82.25913785642317,78.57377420170891,78.0841337171392,77.9169464450095,77.90894306307297,74.47720432824183,73.89586723204289,71.51964410518605,70.32927827695119,70.2258017750965,70.18793173762184,69.99087107813214,69.85942307711323,69.6456277345401,69.55147536261082,69.25393416090424,69.14086375332474,69.0141158548792,67.78378341300639,67.24810324476097,66.78525863074725,66.77975577116197,66.52690337068246,66.12230108383231,66.09131883238062,65.91757385914839,65.7251022318572,65.6433511104564,65.5572881219508,65.20051293168005,65.16378772956423,64.72368724181246,63.78363614146185,62.30449110997799,61.45296728507076,61.01337476048234,60.73368334408425,60.48864151073824,59.065610442186255,58.80309769494216,54.39680124169527,54.25245650872249,53.560738865212485,52.68822141847846,49.49481237650922,48.58382578787084,48.513973787720296,48.30298145651841,46.42957608331301,46.41066289369482,46.3615058994251,46.24403708467408,40.19733654499189,39.393425714136725,39.34772794320706,39.12261331212584,39.08760694159392,38.16488146339577,38.08280966130497,37.949833321222364,34.49296714297766,34.04816268562471,32.77984613475874,30.726208949497632,29.21126778122316,29.022850642484276,28.489974551963993,27.29334572747111,22.394614234240507,21.988532919974254,21.456182591648155,20.980566874180656,20.794801286381936,20.30849345643626,19.74634033109611,19.3880316352086,16.33066901298038,14.109011849235593,13.096867445263703,12.597884421701972,10.878236981995371,7.331909962305469,7.141557023777552,7.123469518316583,7.016493166431487,5.829948710825167,2.8331078875374702,1.5900106632198174,1.542223433042018,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
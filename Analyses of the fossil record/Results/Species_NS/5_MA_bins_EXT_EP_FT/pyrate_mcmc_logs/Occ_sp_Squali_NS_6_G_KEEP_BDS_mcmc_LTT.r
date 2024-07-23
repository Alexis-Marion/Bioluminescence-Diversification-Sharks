

pdf(file='../Results/Species_NS/5_MA_bins_EXT_EP_FT/pyrate_mcmc_logs/Occ_sp_Squali_NS_6_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_6_G_KEEP_BDS_mcmc"
ts=c(46.27905230437605, 65.37879424518027,64.90099607047728,82.53133305698039,85.82569554817047,22.007184203012116,27.197584288265617,6.774960701672513,85.08278568124723,74.60559317310184,6.681641182089028,70.5021086819002,85.19159396685374,23.47128657727972,22.055425893321846,28.6450255608421,54.00075379413578,58.342967405792976,69.2173929917317,84.12366310796202,47.61156457528609,26.159715374714104,72.4259722800028,21.268669225437176,59.3029326435765,64.33987392737333,30.277705316453652,6.791355858068343,21.88085775104654,49.50879215798589,77.77792481575658,71.94001662069645,77.7333180554523,97.87366999395124,109.09881648937612,111.17502757569616,113.42435100170574,90.7509124340269,6.008718703674962,67.91444956040789,49.804668482732524,46.6272667457185,27.73924661424325,75.7264874222492,7.196591836744055,67.00342797809479,33.31461602992355,70.74640360204427,70.1461718390708,64.24733216020483,21.82931904920464,3.8047718947928346,83.54562380270444,53.38696314951614,53.22178185358394,77.9244203205778,22.356220126322345)
te=c(39.37348899420499, 60.64944567230799,63.04714067278551,66.72220908473813,83.10532805446931,0.0,20.118606486282406,0.0,62.412985429604845,65.92604064252707,0.0,65.34592427063404,70.88155616818852,0.0,0.0,18.898447277066623,46.02380423692918,52.078843904265085,66.88533237152691,70.05201855901517,40.61991079762318,15.927876781629676,62.901747606319425,3.0521727449800946,35.39461571483469,49.875385204846026,28.9395195042031,0.0,10.541490006230942,39.28718424673265,65.69360562394786,65.09697265977157,73.96266608351827,93.64730390681682,94.16535320228219,102.78552431271116,94.77057016668952,72.06954956922738,0.0,54.36631468908256,39.18628151945902,38.225738340584044,12.612437544302553,67.9836131867493,0.0,51.30877262459286,12.736119920755657,61.685274557041254,67.25792324064541,39.71856092304357,3.8099299449132538,0.0,66.61969511488665,39.48439513441558,39.34347778005589,72.0964976207183,0.0)
div_traj=c(0, 1,2,3,2,3,2,1,0,1,2,3,4,5,6,5,6,7,8,9,10,11,10,11,10,9,10,9,10,11,12,11,12,11,12,11,12,11,10,9,8,7,8,7,6,7,8,9,8,7,6,5,4,5,6,5,6,7,8,7,6,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1,0,1,2,1,2,3,4,5,6,7,8,9,10,11,12,11,10,9,8,7,6,7,8,9,10,11,10,11,10)
time_events=c(113.42435100170574, 111.17502757569616,109.09881648937612,102.78552431271116,97.87366999395124,94.77057016668952,94.16535320228219,93.64730390681682,90.7509124340269,85.82569554817047,85.19159396685374,85.08278568124723,84.12366310796202,83.54562380270444,83.10532805446931,82.53133305698039,77.9244203205778,77.77792481575658,77.7333180554523,75.7264874222492,74.60559317310184,73.96266608351827,72.4259722800028,72.0964976207183,72.06954956922738,71.94001662069645,70.88155616818852,70.74640360204427,70.5021086819002,70.1461718390708,70.05201855901517,69.2173929917317,67.9836131867493,67.91444956040789,67.25792324064541,67.00342797809479,66.88533237152691,66.72220908473813,66.61969511488665,65.92604064252707,65.69360562394786,65.37879424518027,65.34592427063404,65.09697265977157,64.90099607047728,64.33987392737333,64.24733216020483,63.04714067278551,62.901747606319425,62.412985429604845,61.685274557041254,60.64944567230799,59.3029326435765,58.342967405792976,54.36631468908256,54.00075379413578,53.38696314951614,53.22178185358394,52.078843904265085,51.30877262459286,49.875385204846026,49.804668482732524,49.50879215798589,47.61156457528609,46.6272667457185,46.27905230437605,46.02380423692918,40.61991079762318,39.71856092304357,39.48439513441558,39.37348899420499,39.34347778005589,39.28718424673265,39.18628151945902,38.225738340584044,35.39461571483469,33.31461602992355,30.277705316453652,28.9395195042031,28.6450255608421,27.73924661424325,27.197584288265617,26.159715374714104,23.47128657727972,22.356220126322345,22.055425893321846,22.007184203012116,21.88085775104654,21.82931904920464,21.268669225437176,20.118606486282406,18.898447277066623,15.927876781629676,12.736119920755657,12.612437544302553,10.541490006230942,7.196591836744055,6.791355858068343,6.774960701672513,6.681641182089028,6.008718703674962,3.8099299449132538,3.8047718947928346,3.0521727449800946,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
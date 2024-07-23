

pdf(file='../Results/Species/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_sp_Squali_17_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_17_G_KEEP_BDS_mcmc"
ts=c(45.06006035945063, 45.62188431881401,13.744883138925687,65.79048256524773,66.56586354655792,84.2376290432852,84.92623090377275,19.331000003881382,18.926227663047886,25.55029405526717,6.306705115187307,29.03148432803419,84.13420360866776,77.76532256046946,45.289664977014155,6.541475700814885,77.58175624526427,85.99236081267547,85.42023893503291,69.64933209547544,67.93191633949915,82.85796282940422,5.8913397298980215,20.92866501206842,15.143831500143985,59.03649729664535,45.371848841496345,22.073720070622095,18.503173019405978,29.293965150326205,54.12900924505914,57.36016631494394,69.82401896312156,84.18345986152897,45.62705898057714,19.46361685323231,29.661010096107155,45.27168457895581,16.806373960079263,75.12082254605906,69.40429306763097,83.9156199597147,49.15245053653508,20.553535218703235,56.85530420973985,64.44294986164365,70.77468143323352,32.159890355680105,6.821172427276553,21.633615549148182,22.366176355968125,49.477637490869654,77.25039865370637,69.21975818944934,94.95858652459107,74.8684370277769,127.04033970126069,100.41713124113964,108.11126072804636,110.26887533389066,111.5644713414751,66.80960995779598,91.93083673490436,77.282072041237,6.059483607645414,24.177992655048694,45.43491493956781,6.610366838649739,17.245895197784304,6.252631759414255,30.232704474792204,6.318304377810702,65.5047862560278,76.25041474473663,48.737198792202136,45.689769465490365,33.66225743488682,6.159278107734947,15.020059553862882,67.31385054368185,78.12610903775312,77.0050697536218,6.143887265739753,65.75407999740085,34.17980676583795,67.21256539228938,69.9423484738419,64.57610468872561,73.84008518304157,21.819944648146098,52.00133344481066,3.4152704309826025,66.37423051008204,82.00932814976325,54.27358242814804,54.2762921492603,75.96044461096652,45.227724995749526,19.606140560007272)
te=c(41.636587172753, 37.74716834031803,11.56126124094239,61.571743361619,64.06898509088826,67.23701759165586,83.11426356389195,0.0,15.522168061549602,14.139546382649568,0.0,26.651247322607087,60.92393848248387,65.46940092570777,41.6152412393407,0.0,74.60406040715682,84.69196277470013,84.06940249903059,66.1212390098693,66.97473424571432,72.66586188358636,2.9220575498979082,0.0,12.803115019302387,57.58955582638183,42.00669931120295,0.0,14.928199408179253,19.485365909534796,48.705736986397774,49.0514766377576,67.25380523165745,69.48343638545377,39.89018469886356,16.061486430063017,18.295870441951692,41.87717201724588,14.16120160312773,72.5264291866386,64.78635688719781,82.47547340243581,46.11343335225119,2.4379186738850853,36.30438787396998,52.48904332041828,69.46979798807635,27.097619252738816,0.0,19.622881189124186,12.389302736659273,37.66653686917677,66.13865908963285,65.55457218356023,94.26736076100477,71.06427584242158,125.8275036146708,93.36679099627588,93.66202789430787,99.80957888712258,94.69595981105441,65.90014779661571,72.24472981934062,74.4992174999031,0.0,22.460682348967467,41.86729591673766,0.0,14.3432297409741,0.0,27.953061198074415,0.0,60.66646217512823,74.0177736915524,37.54008311626752,39.58119157538693,12.948257642238088,0.0,12.592757435203058,66.41634975135774,69.10667808863262,74.31439816668123,0.0,48.113265527538665,13.612171754484562,62.60166171865481,66.99605521276402,41.88756833889909,71.349161760357,1.9648610456654505,48.74545468200823,0.0,65.38279089400112,69.2837068874462,40.171647723563076,39.620150337730735,72.86996723362236,41.86737904373053,0.0)
div_traj=c(0, 1,0,1,2,3,4,3,4,3,2,1,0,1,2,3,4,3,4,5,6,5,6,5,6,5,6,7,8,9,10,11,12,13,14,15,16,15,14,13,12,13,12,11,10,9,8,7,8,9,10,11,10,9,10,9,10,9,10,11,10,9,10,9,8,9,10,9,10,9,8,7,8,9,8,9,8,7,6,7,8,7,6,5,4,3,4,3,4,5,6,7,8,7,8,9,10,9,8,9,8,7,6,7,8,9,10,11,12,13,14,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,6,5,4,5,6,5,6,7,8,9,10,11,10,11,10,11,12,13,14,13,14,15,14,13,14,15,14,13,12,11,12,11,10,9,8,7,6,7,8,9,10,11,12,13,14,15,16,17,16,15,14)
time_events=c(127.04033970126069, 125.8275036146708,111.5644713414751,110.26887533389066,108.11126072804636,100.41713124113964,99.80957888712258,94.95858652459107,94.69595981105441,94.26736076100477,93.66202789430787,93.36679099627588,91.93083673490436,85.99236081267547,85.42023893503291,84.92623090377275,84.69196277470013,84.2376290432852,84.18345986152897,84.13420360866776,84.06940249903059,83.9156199597147,83.11426356389195,82.85796282940422,82.47547340243581,82.00932814976325,78.12610903775312,77.76532256046946,77.58175624526427,77.282072041237,77.25039865370637,77.0050697536218,76.25041474473663,75.96044461096652,75.12082254605906,74.8684370277769,74.60406040715682,74.4992174999031,74.31439816668123,74.0177736915524,73.84008518304157,72.86996723362236,72.66586188358636,72.5264291866386,72.24472981934062,71.349161760357,71.06427584242158,70.77468143323352,69.9423484738419,69.82401896312156,69.64933209547544,69.48343638545377,69.46979798807635,69.40429306763097,69.2837068874462,69.21975818944934,69.10667808863262,67.93191633949915,67.31385054368185,67.25380523165745,67.23701759165586,67.21256539228938,66.99605521276402,66.97473424571432,66.80960995779598,66.56586354655792,66.41634975135774,66.37423051008204,66.13865908963285,66.1212390098693,65.90014779661571,65.79048256524773,65.75407999740085,65.55457218356023,65.5047862560278,65.46940092570777,65.38279089400112,64.78635688719781,64.57610468872561,64.44294986164365,64.06898509088826,62.60166171865481,61.571743361619,60.92393848248387,60.66646217512823,59.03649729664535,57.58955582638183,57.36016631494394,56.85530420973985,54.2762921492603,54.27358242814804,54.12900924505914,52.48904332041828,52.00133344481066,49.477637490869654,49.15245053653508,49.0514766377576,48.74545468200823,48.737198792202136,48.705736986397774,48.113265527538665,46.11343335225119,45.689769465490365,45.62705898057714,45.62188431881401,45.43491493956781,45.371848841496345,45.289664977014155,45.27168457895581,45.227724995749526,45.06006035945063,42.00669931120295,41.88756833889909,41.87717201724588,41.86737904373053,41.86729591673766,41.636587172753,41.6152412393407,40.171647723563076,39.89018469886356,39.620150337730735,39.58119157538693,37.74716834031803,37.66653686917677,37.54008311626752,36.30438787396998,34.17980676583795,33.66225743488682,32.159890355680105,30.232704474792204,29.661010096107155,29.293965150326205,29.03148432803419,27.953061198074415,27.097619252738816,26.651247322607087,25.55029405526717,24.177992655048694,22.460682348967467,22.366176355968125,22.073720070622095,21.819944648146098,21.633615549148182,20.92866501206842,20.553535218703235,19.622881189124186,19.606140560007272,19.485365909534796,19.46361685323231,19.331000003881382,18.926227663047886,18.503173019405978,18.295870441951692,17.245895197784304,16.806373960079263,16.061486430063017,15.522168061549602,15.143831500143985,15.020059553862882,14.928199408179253,14.3432297409741,14.16120160312773,14.139546382649568,13.744883138925687,13.612171754484562,12.948257642238088,12.803115019302387,12.592757435203058,12.389302736659273,11.56126124094239,6.821172427276553,6.610366838649739,6.541475700814885,6.318304377810702,6.306705115187307,6.252631759414255,6.159278107734947,6.143887265739753,6.059483607645414,5.8913397298980215,3.4152704309826025,2.9220575498979082,2.4379186738850853,1.9648610456654505,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                


pdf(file='../Results/Species_NS/5_MA_bins_EP/pyrate_mcmc_logs/Occ_sp_Squali_NS_20_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_20_G_KEEP_BDS_mcmc"
ts=c(46.535428907274905, 65.70764728263899,67.21729120725207,84.18978771049476,85.4950394936273,19.20969552349867,28.543321442733195,5.761428782536738,84.81329127117179,76.72682626469575,6.041998697104908,71.0613466751282,85.10981338739897,22.557217748160557,21.49293345433136,26.526422470873335,53.31466981292006,58.02272663833031,69.71236254471,83.91046897547675,46.06396283044487,33.82877995843006,69.21925226222649,21.056682331590768,58.45931111644581,66.14276554171562,33.92104851012806,6.795874290188002,22.659902043381997,48.696411511171476,78.34070196620397,68.25753554044212,78.14113642082793,98.50305776602308,108.45388556987415,109.42052005910041,111.84986801744024,91.25748267894576,5.628253174310872,66.4410174786288,49.098809517354425,46.294257639492095,30.006542957313684,76.95541085430757,6.513826270717095,66.56256902716514,33.997309511646485,69.52605386582901,70.83881353881743,66.58486908193557,18.339488642131936,3.3404787042864976,84.96106415346712,54.23612619104271,53.005883824746824,75.56872032759905,20.43681423499441)
te=c(37.741981777026076, 60.91225685174129,62.82316752383892,65.55090272807861,83.59119988644547,0.0,16.121508714943424,0.0,62.53759595142296,66.20403269269863,0.0,65.83550788744157,71.96228462487535,0.0,0.0,21.235200045302044,49.270425725764426,49.79769894946074,65.62035276304424,69.24710569916347,40.584131408029684,13.5245708042775,61.08530270959577,2.6094785474932474,33.74836344957666,44.77228196501109,29.292178980204916,0.0,11.940381503235105,39.262966709835965,65.33644200988338,66.66864079413968,68.19957886099431,94.1683592772471,93.76258435916493,102.81821257862609,94.96408488545316,73.16753710528921,0.0,53.10109856659473,39.581111960578845,38.30557810468206,10.377017264399706,66.03469735489175,0.0,53.321762365082066,13.065690906267404,61.93052543379495,67.61040110337662,39.83280905238042,3.292073429145538,0.0,69.8487985042853,39.52932692378536,39.11265525872428,71.78209239914048,0.0)
div_traj=c(0, 1,2,3,2,3,2,1,0,1,2,3,4,5,6,7,6,7,8,9,10,11,10,9,8,9,10,9,10,11,10,11,12,11,10,11,10,11,12,13,12,13,12,11,12,11,10,9,8,7,6,5,4,5,6,7,6,7,6,7,6,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1,2,3,4,3,4,3,4,5,6,7,8,7,8,9,10,11,10,9,8,7,6,7,8,9,10,11,12,11,10)
time_events=c(111.84986801744024, 109.42052005910041,108.45388556987415,102.81821257862609,98.50305776602308,94.96408488545316,94.1683592772471,93.76258435916493,91.25748267894576,85.4950394936273,85.10981338739897,84.96106415346712,84.81329127117179,84.18978771049476,83.91046897547675,83.59119988644547,78.34070196620397,78.14113642082793,76.95541085430757,76.72682626469575,75.56872032759905,73.16753710528921,71.96228462487535,71.78209239914048,71.0613466751282,70.83881353881743,69.8487985042853,69.71236254471,69.52605386582901,69.24710569916347,69.21925226222649,68.25753554044212,68.19957886099431,67.61040110337662,67.21729120725207,66.66864079413968,66.58486908193557,66.56256902716514,66.4410174786288,66.20403269269863,66.14276554171562,66.03469735489175,65.83550788744157,65.70764728263899,65.62035276304424,65.55090272807861,65.33644200988338,62.82316752383892,62.53759595142296,61.93052543379495,61.08530270959577,60.91225685174129,58.45931111644581,58.02272663833031,54.23612619104271,53.321762365082066,53.31466981292006,53.10109856659473,53.005883824746824,49.79769894946074,49.270425725764426,49.098809517354425,48.696411511171476,46.535428907274905,46.294257639492095,46.06396283044487,44.77228196501109,40.584131408029684,39.83280905238042,39.581111960578845,39.52932692378536,39.262966709835965,39.11265525872428,38.30557810468206,37.741981777026076,33.997309511646485,33.92104851012806,33.82877995843006,33.74836344957666,30.006542957313684,29.292178980204916,28.543321442733195,26.526422470873335,22.659902043381997,22.557217748160557,21.49293345433136,21.235200045302044,21.056682331590768,20.43681423499441,19.20969552349867,18.339488642131936,16.121508714943424,13.5245708042775,13.065690906267404,11.940381503235105,10.377017264399706,6.795874290188002,6.513826270717095,6.041998697104908,5.761428782536738,5.628253174310872,3.3404787042864976,3.292073429145538,2.6094785474932474,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
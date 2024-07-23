

pdf(file='../Results/Species_NS/5_MA_bins_EXT_EP_FT/pyrate_mcmc_logs/Occ_sp_Squali_NS_20_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_20_G_KEEP_BDS_mcmc"
ts=c(46.70958277500198, 65.82471816991809,67.3243875827018,84.28361693886917,85.55169151224506,19.695121610269528,28.81301798614754,5.8307574553186505,84.74554467690189,76.98072551196292,5.978486764496752,71.09259611662588,85.08726240520662,22.756809087354185,21.858555474543497,26.4953742184327,53.43899230907763,58.13021383689444,69.92207541399723,83.97668170338801,46.118652762079684,33.509556902488924,69.26116814800072,21.121181184541477,58.48637752844127,66.73172742170898,33.735003912559215,7.300617876834857,22.81252645842292,48.74456455615495,78.5418619312293,68.38699008497792,78.46806633465047,98.35240276627088,108.4010424739223,109.324605373577,111.7431864375224,91.25637591398566,5.369527410932689,67.37681677469209,49.112175582135976,46.511379302219396,30.144116812818403,77.0371196753399,7.0140888868218125,66.90614981679552,33.80564860088678,69.60061916221042,70.89459319060502,67.002657913034,19.097793626986732,2.4443241030448273,85.02150269761994,54.31633929096229,53.02365701957796,75.75571716538032,20.86338026771048)
te=c(37.687307325568774, 60.93491715793089,62.84341662129099,65.52422675702962,83.5802861062985,0.0,15.69548127024173,0.0,62.67149682152953,66.07505296693736,0.0,65.75892940055752,71.8930745724651,0.0,0.0,21.20184576426615,48.948571530099464,49.86074994273397,65.49302852934794,69.26672673212772,40.560690307654546,13.390482460082024,61.12343793231566,2.789457243743294,34.08501380190182,44.477906002720395,29.110664434686505,0.0,11.920434564390655,39.273196041136714,65.25232784318808,66.5613261849071,68.09039407016235,94.16305015774216,93.73431859916734,102.69458768720358,94.89752931698041,73.03012548502448,0.0,53.251034982428365,39.47351631090976,38.25690268287159,10.501019037512998,65.96650903673827,0.0,53.15171910618159,13.020529400686073,62.00976991140656,67.45064427019997,39.872531391848874,3.6444211424154083,0.0,69.84243796971285,39.4460632750294,39.16791539805219,71.61548639336377,0.0)
div_traj=c(0, 1,2,3,2,3,2,1,0,1,2,3,4,5,6,7,6,7,8,9,10,11,10,9,8,9,10,11,10,11,10,11,12,11,10,11,12,13,14,15,14,13,12,13,12,11,10,9,8,7,6,5,4,5,6,7,8,7,6,7,6,7,6,7,8,9,10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,3,4,5,6,7,8,7,8,9,10,11,10,9,8,7,6,7,8,9,10,11,10,9,10)
time_events=c(111.7431864375224, 109.324605373577,108.4010424739223,102.69458768720358,98.35240276627088,94.89752931698041,94.16305015774216,93.73431859916734,91.25637591398566,85.55169151224506,85.08726240520662,85.02150269761994,84.74554467690189,84.28361693886917,83.97668170338801,83.5802861062985,78.5418619312293,78.46806633465047,77.0371196753399,76.98072551196292,75.75571716538032,73.03012548502448,71.8930745724651,71.61548639336377,71.09259611662588,70.89459319060502,69.92207541399723,69.84243796971285,69.60061916221042,69.26672673212772,69.26116814800072,68.38699008497792,68.09039407016235,67.45064427019997,67.37681677469209,67.3243875827018,67.002657913034,66.90614981679552,66.73172742170898,66.5613261849071,66.07505296693736,65.96650903673827,65.82471816991809,65.75892940055752,65.52422675702962,65.49302852934794,65.25232784318808,62.84341662129099,62.67149682152953,62.00976991140656,61.12343793231566,60.93491715793089,58.48637752844127,58.13021383689444,54.31633929096229,53.43899230907763,53.251034982428365,53.15171910618159,53.02365701957796,49.86074994273397,49.112175582135976,48.948571530099464,48.74456455615495,46.70958277500198,46.511379302219396,46.118652762079684,44.477906002720395,40.560690307654546,39.872531391848874,39.47351631090976,39.4460632750294,39.273196041136714,39.16791539805219,38.25690268287159,37.687307325568774,34.08501380190182,33.80564860088678,33.735003912559215,33.509556902488924,30.144116812818403,29.110664434686505,28.81301798614754,26.4953742184327,22.81252645842292,22.756809087354185,21.858555474543497,21.20184576426615,21.121181184541477,20.86338026771048,19.695121610269528,19.097793626986732,15.69548127024173,13.390482460082024,13.020529400686073,11.920434564390655,10.501019037512998,7.300617876834857,7.0140888868218125,5.978486764496752,5.8307574553186505,5.369527410932689,3.6444211424154083,2.789457243743294,2.4443241030448273,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
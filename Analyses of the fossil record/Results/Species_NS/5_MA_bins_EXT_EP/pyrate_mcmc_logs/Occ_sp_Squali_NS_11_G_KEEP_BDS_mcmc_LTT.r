

pdf(file='../Results/Species_NS/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_sp_Squali_NS_11_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_11_G_KEEP_BDS_mcmc"
ts=c(46.03568608040039, 63.76131908763809,66.96401930993888,84.32145358378325,86.43474064444312,19.271701660279692,29.642548108049684,6.499388011612549,85.24154822225957,76.9835717297168,6.103760797064394,71.14755133080344,85.15249508858804,22.564705031213283,21.54715196911115,29.56928599708604,53.752222148006744,58.287039043921155,70.98130633365207,84.34553146702834,46.541910530551164,33.60413367936,67.87420132297358,21.153452635410883,57.99728272620176,62.27699178184914,33.88595943395474,6.701588823212226,21.805215796064353,48.431615721726956,76.0964009156741,70.55106428309358,84.78203517781395,99.97609007377649,104.41669226617337,111.18959504637158,111.43677345938488,91.32043164951534,6.22345728087018,67.54392839039942,48.33503462009518,46.26830172270689,30.584590140040593,75.49312012228164,6.212112737699267,67.14296685162986,33.789449269784505,71.28287540024124,70.22197218761103,66.64702970152715,21.171843766041896,3.9984408956075193,83.86562215211461,54.38336113164012,54.237453195795354,78.81386709383705,21.202111223126494)
te=c(39.493504092025596, 62.16616774568509,60.64123344303687,66.1304953439515,81.86534361810759,0.0,19.734331460110123,0.0,62.000539242801516,65.473760083342,0.0,66.50519938789232,72.07627612362987,0.0,0.0,17.5214762385414,48.06970814274491,52.03730591489612,65.82183542751942,69.63870502372833,40.31217127299315,13.50456484340422,65.06800399041428,2.5884551594473155,37.32561811452333,49.47728866483942,30.75230169514568,0.0,12.380305986729768,39.46592998664072,65.14675922808951,65.03537220248117,73.10072994342109,94.55565232210144,94.55105750257096,105.31409695192666,94.80303739287086,72.86927699385176,0.0,53.01460096229184,39.188955470606516,38.79614383562243,12.215632894754778,68.09663314950011,0.0,51.17227155375036,12.675087800269782,61.051575290524575,66.35628714035242,39.52319584987724,3.1979108511825425,0.0,69.07820997074498,39.45339920280851,39.74992541718967,70.99823574315978,0.0)
div_traj=c(0, 1,2,1,2,3,2,1,0,1,2,3,4,5,6,7,8,7,8,9,10,11,10,9,8,9,10,9,10,11,12,11,10,9,10,11,12,13,14,13,12,11,10,9,8,7,6,7,8,7,6,5,4,5,6,7,8,9,8,7,6,5,6,7,6,7,8,9,8,7,6,5,4,3,2,1,0,1,2,3,2,3,4,5,6,7,8,9,10,11,10,11,10,9,8,7,6,7,8,9,10,11,12,11,10)
time_events=c(111.43677345938488, 111.18959504637158,105.31409695192666,104.41669226617337,99.97609007377649,94.80303739287086,94.55565232210144,94.55105750257096,91.32043164951534,86.43474064444312,85.24154822225957,85.15249508858804,84.78203517781395,84.34553146702834,84.32145358378325,83.86562215211461,81.86534361810759,78.81386709383705,76.9835717297168,76.0964009156741,75.49312012228164,73.10072994342109,72.86927699385176,72.07627612362987,71.28287540024124,71.14755133080344,70.99823574315978,70.98130633365207,70.55106428309358,70.22197218761103,69.63870502372833,69.07820997074498,68.09663314950011,67.87420132297358,67.54392839039942,67.14296685162986,66.96401930993888,66.64702970152715,66.50519938789232,66.35628714035242,66.1304953439515,65.82183542751942,65.473760083342,65.14675922808951,65.06800399041428,65.03537220248117,63.76131908763809,62.27699178184914,62.16616774568509,62.000539242801516,61.051575290524575,60.64123344303687,58.287039043921155,57.99728272620176,54.38336113164012,54.237453195795354,53.752222148006744,53.01460096229184,52.03730591489612,51.17227155375036,49.47728866483942,48.431615721726956,48.33503462009518,48.06970814274491,46.541910530551164,46.26830172270689,46.03568608040039,40.31217127299315,39.74992541718967,39.52319584987724,39.493504092025596,39.46592998664072,39.45339920280851,39.188955470606516,38.79614383562243,37.32561811452333,33.88595943395474,33.789449269784505,33.60413367936,30.75230169514568,30.584590140040593,29.642548108049684,29.56928599708604,22.564705031213283,21.805215796064353,21.54715196911115,21.202111223126494,21.171843766041896,21.153452635410883,19.734331460110123,19.271701660279692,17.5214762385414,13.50456484340422,12.675087800269782,12.380305986729768,12.215632894754778,6.701588823212226,6.499388011612549,6.22345728087018,6.212112737699267,6.103760797064394,3.9984408956075193,3.1979108511825425,2.5884551594473155,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
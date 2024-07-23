

pdf(file='../Results/Species/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_sp_Squali_9_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_9_G_KEEP_BDS_mcmc"
ts=c(44.24559962975814, 45.27384313183361,16.409192639033982,65.67137322264064,66.7539099882522,82.73504138798477,85.38591763781983,18.485197966505783,20.71016106819047,28.539512842742944,6.165377541731531,26.006740968359185,83.2400266678149,77.48524586111458,44.27970270379387,6.44897840366416,74.54800910384498,86.02207007540039,84.9156795981643,70.28945305230397,83.33723048675306,83.06351343318806,5.717170003856599,23.47296373424685,13.719429164698036,58.70294671297735,44.222534637044035,17.916255187283234,17.171258862644958,28.609213565334628,53.0321865960182,60.02890516606909,70.40022343572443,84.25694002475501,44.76002430324576,16.902529721865587,27.98063090117153,44.18950994238835,16.597528607279305,74.90085925179419,69.01057354197803,82.35274955804158,62.78483153201541,20.97868211234414,58.86227452916558,66.94389614808196,72.40452022759192,32.39869995021397,6.387109939157246,22.253196448017512,24.280081665698653,50.04540842358352,76.71594229009179,69.52736251433157,94.72433203776947,75.3377969363496,127.65064745454762,101.32066195249006,103.77723382718837,109.6793360627438,113.9182501949752,70.10934844881186,90.36015645227273,73.80728073427095,6.655990324112354,27.70355519711126,44.06755237223988,5.869429459690062,16.73384906249241,6.181929967004732,29.99312234269409,6.3103755527416645,65.84904293914661,73.91820430936095,50.03864070126142,44.50891617779445,25.349471884653862,6.248187265424978,14.405490197001392,70.82184789575135,76.23240271470593,74.63686413104587,6.561397236438381,66.8090992195392,34.44862401798432,70.16891308078918,72.01687092919143,65.03385710293158,77.32308958974095,12.75206132926942,53.020288007740206,2.688200015929849,66.298707254595,82.2244712280072,53.864173063539965,53.412149363471826,77.50092115725992,44.45579457533313,19.849968250032944)
te=c(42.25682086128319, 40.28769140195611,14.438416042205276,61.299483728531456,63.508486790798806,65.89901170050807,83.46100606220715,0.0,18.952391066931703,21.08393467377379,0.0,24.74863822160765,61.63777600848551,65.8911153911894,42.28337319680631,0.0,72.67683270314318,84.88645801397975,83.6125231927732,65.69786350463549,82.11443109934899,71.69339091520513,2.8318848207702185,0.0,11.706808977140966,56.29037165556868,42.321291412430185,0.0,15.026025769111676,21.20290302775967,48.32186714068084,46.175109505352374,66.75538877662069,70.26290490338668,41.310755217757475,14.88614061598423,21.180669467698074,42.14981890222382,14.515561890243362,73.28544981016313,64.96718162640627,81.15805926020545,61.53853101222042,3.3727146187439834,36.75384498464518,49.62253036464591,70.79491803369109,29.148850583152683,0.0,20.52008122428664,11.169038229664846,40.37770797459127,67.42154752737014,65.68272738850926,94.48170433067331,65.72808144453636,126.41242862630008,94.80213328700073,94.5697134889362,104.70384268922116,94.94530579543283,68.33984579070406,73.99058081884786,72.1221907947716,0.0,26.388911369448103,42.1136869411976,0.0,14.626615715603698,0.0,28.719866424564184,0.0,54.7906565640664,72.17713373663422,38.55546490817017,40.054926242100116,12.772085860916057,0.0,12.259818623353922,69.31826211359075,65.73394729959432,72.90695233126772,0.0,47.80375793485461,13.145347671592509,61.271461189205645,68.20429991432437,41.056608141589834,74.283658322874,2.897860649994617,49.99687610843678,0.0,64.75336885580921,70.55217452373178,40.27618572461343,40.82919501391396,70.76874514675818,42.43174794966537,0.0)
div_traj=c(0, 1,0,1,2,1,2,3,2,1,2,1,0,1,2,3,4,3,4,3,2,3,4,5,6,7,8,7,6,7,8,9,10,11,12,13,14,15,14,13,14,15,14,13,12,13,12,11,12,11,12,11,10,9,10,11,10,11,12,13,12,13,12,11,10,11,12,11,12,13,12,11,12,11,10,9,8,9,10,9,8,7,8,7,6,5,4,5,6,7,6,5,6,7,8,9,10,11,10,9,8,7,6,7,8,9,10,11,12,13,14,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,1,2,3,2,1,2,3,4,5,4,5,6,5,6,7,8,7,6,5,6,7,6,7,6,7,8,9,10,11,12,13,12,11,10,9,8,9,10,9,8,9,8,7,6,7,8,9,10,11,12,13,14,15,16,15,14,13,14)
time_events=c(127.65064745454762, 126.41242862630008,113.9182501949752,109.6793360627438,104.70384268922116,103.77723382718837,101.32066195249006,94.94530579543283,94.80213328700073,94.72433203776947,94.5697134889362,94.48170433067331,90.36015645227273,86.02207007540039,85.38591763781983,84.9156795981643,84.88645801397975,84.25694002475501,83.6125231927732,83.46100606220715,83.33723048675306,83.2400266678149,83.06351343318806,82.73504138798477,82.35274955804158,82.2244712280072,82.11443109934899,81.15805926020545,77.50092115725992,77.48524586111458,77.32308958974095,76.71594229009179,76.23240271470593,75.3377969363496,74.90085925179419,74.63686413104587,74.54800910384498,74.283658322874,73.99058081884786,73.91820430936095,73.80728073427095,73.28544981016313,72.90695233126772,72.67683270314318,72.40452022759192,72.17713373663422,72.1221907947716,72.01687092919143,71.69339091520513,70.82184789575135,70.79491803369109,70.76874514675818,70.55217452373178,70.40022343572443,70.28945305230397,70.26290490338668,70.16891308078918,70.10934844881186,69.52736251433157,69.31826211359075,69.01057354197803,68.33984579070406,68.20429991432437,67.42154752737014,66.94389614808196,66.8090992195392,66.75538877662069,66.7539099882522,66.298707254595,65.89901170050807,65.8911153911894,65.84904293914661,65.73394729959432,65.72808144453636,65.69786350463549,65.68272738850926,65.67137322264064,65.03385710293158,64.96718162640627,64.75336885580921,63.508486790798806,62.78483153201541,61.63777600848551,61.53853101222042,61.299483728531456,61.271461189205645,60.02890516606909,58.86227452916558,58.70294671297735,56.29037165556868,54.7906565640664,53.864173063539965,53.412149363471826,53.0321865960182,53.020288007740206,50.04540842358352,50.03864070126142,49.99687610843678,49.62253036464591,48.32186714068084,47.80375793485461,46.175109505352374,45.27384313183361,44.76002430324576,44.50891617779445,44.45579457533313,44.27970270379387,44.24559962975814,44.222534637044035,44.18950994238835,44.06755237223988,42.43174794966537,42.321291412430185,42.28337319680631,42.25682086128319,42.14981890222382,42.1136869411976,41.310755217757475,41.056608141589834,40.82919501391396,40.37770797459127,40.28769140195611,40.27618572461343,40.054926242100116,38.55546490817017,36.75384498464518,34.44862401798432,32.39869995021397,29.99312234269409,29.148850583152683,28.719866424564184,28.609213565334628,28.539512842742944,27.98063090117153,27.70355519711126,26.388911369448103,26.006740968359185,25.349471884653862,24.74863822160765,24.280081665698653,23.47296373424685,22.253196448017512,21.20290302775967,21.180669467698074,21.08393467377379,20.97868211234414,20.71016106819047,20.52008122428664,19.849968250032944,18.952391066931703,18.485197966505783,17.916255187283234,17.171258862644958,16.902529721865587,16.73384906249241,16.597528607279305,16.409192639033982,15.026025769111676,14.88614061598423,14.626615715603698,14.515561890243362,14.438416042205276,14.405490197001392,13.719429164698036,13.145347671592509,12.772085860916057,12.75206132926942,12.259818623353922,11.706808977140966,11.169038229664846,6.655990324112354,6.561397236438381,6.44897840366416,6.387109939157246,6.3103755527416645,6.248187265424978,6.181929967004732,6.165377541731531,5.869429459690062,5.717170003856599,3.3727146187439834,2.897860649994617,2.8318848207702185,2.688200015929849,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
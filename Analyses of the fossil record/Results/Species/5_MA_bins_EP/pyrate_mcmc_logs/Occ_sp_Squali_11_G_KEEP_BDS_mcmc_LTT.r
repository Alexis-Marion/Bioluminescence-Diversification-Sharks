

pdf(file='../Results/Species/5_MA_bins_EP/pyrate_mcmc_logs/Occ_sp_Squali_11_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_11_G_KEEP_BDS_mcmc"
ts=c(44.729707097891634, 45.47670698823559,13.357781713932557,65.18833104589127,64.01338306567423,83.82774277670238,86.15761358910542,21.36027472173769,20.774202827642693,25.38273514075396,6.41422081661375,24.67652712609007,84.86694214381284,74.58361820033605,45.211042509769136,6.323430768076129,77.3858228366621,85.94649168983183,84.92853925369353,70.59489107115465,83.88171428804769,83.94730150043614,5.988978382332227,20.052440768387406,12.936051624355805,59.501725048668156,44.77623598436978,19.35892677385418,20.852247937344746,25.490008250305177,52.70474076045028,58.742413649236255,70.80725840033709,82.05726016290478,45.41839217187456,18.35686013324457,27.885322477001658,44.62659634607659,17.29334843585728,74.81604018387954,67.57031378183422,84.16486234950084,61.68103571158919,21.09349463960433,59.60258979270896,65.2060580294703,69.53117718125154,33.24546949495093,6.566564032023734,23.577190034623282,23.83479381888901,49.40084938640584,75.84620494320025,71.28027873049473,95.14122549009053,74.35670041501308,129.01500585432686,99.37079484077157,105.91301845831713,109.39277729695725,111.887157704515,70.10402260920615,91.22987438771702,73.11744465777583,5.948105629162042,26.89767806161609,45.29854611720612,5.953663700841609,17.27460899582585,6.18927334696072,32.04444839294616,6.087052810453401,65.81837592333682,76.57026961039973,48.38536052890897,45.15247606786069,27.010938977679114,6.2725803662237425,12.693293877613826,67.84681011403126,74.56843162508834,76.12284670059285,6.540356639138789,66.22201681367976,33.65516414161626,70.8036729100543,67.49210331323822,66.71306083916912,75.04984429977227,20.957866383366827,54.32731105384171,4.899697738441464,62.96151654442945,84.48672796708584,53.469482771859724,53.970529417556335,75.32109008676477,44.938658094208186,23.026230758386063)
te=c(42.10547886717194, 37.22262153464601,11.5622321506254,63.55208562562408,62.15693255499719,65.72072397937451,83.55975235164334,0.0,18.092108649551577,15.461985061347962,0.0,22.89963958324237,61.12846369624528,64.91381911031287,42.50001401067301,0.0,74.51867261711911,84.62791851948175,83.41920860330092,66.93047730912309,82.36944496088982,71.97577532822777,2.8180585211394806,0.0,11.252575124466176,57.49394383110114,42.245345030701486,0.0,18.6977980136776,21.556648568603347,48.442452387338925,45.84403019292141,68.7583591856931,70.32924850620888,42.25147248296379,15.456634123225468,16.119283038362738,42.02962522632461,14.738922931615443,73.12545115215933,63.408774691407906,82.59985442166085,60.53446681780395,2.2537574560875893,34.051162682640495,45.41866897559957,68.11043541442514,29.798721826132354,0.0,21.725530299425216,12.085153010540349,38.994746522073655,65.99828309402207,66.19539050842751,94.46694310578393,65.45091124632094,127.73636670322533,93.73874634479658,93.43505788773689,102.74257080475182,94.77385445877509,68.76545788465243,73.05191117385225,71.49005471145058,0.0,25.570356360313852,42.56847210635154,0.0,14.619956066410086,0.0,31.15093598738311,0.0,54.29661161398992,74.14454483801123,40.05187809911719,39.22623428954331,12.295631198551366,0.0,11.241104728298867,66.39571825493401,70.5137912923328,74.19573987084722,0.0,50.140369305066365,14.213845410890979,60.94040505871836,65.90381340903588,41.292287162752835,73.44516234286183,3.344716682179103,51.38101859508939,0.0,61.704988716942594,69.87603030464756,39.44116550297453,40.48369623343603,71.37571241646155,42.38038561219,0.0)
div_traj=c(0, 1,0,1,2,3,2,3,4,3,2,1,0,1,2,3,4,5,4,5,6,7,8,9,8,7,6,5,6,7,8,9,10,11,12,13,14,15,14,15,14,13,12,11,12,11,10,9,8,9,10,11,12,11,10,11,10,11,10,9,8,9,10,11,10,11,10,11,10,9,8,9,8,7,8,9,8,9,8,7,8,7,6,7,6,5,4,5,6,7,6,7,6,7,8,9,8,7,8,7,8,7,8,7,8,9,10,11,12,13,14,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,1,2,3,2,1,2,3,4,3,4,5,6,7,8,9,8,7,6,7,8,9,10,11,12,13,12,13,12,13,14,13,12,11,10,9,8,9,10,11,10,9,8,7,6,7,8,9,10,11,12,13,14,15,16,17,16,15,14)
time_events=c(129.01500585432686, 127.73636670322533,111.887157704515,109.39277729695725,105.91301845831713,102.74257080475182,99.37079484077157,95.14122549009053,94.77385445877509,94.46694310578393,93.73874634479658,93.43505788773689,91.22987438771702,86.15761358910542,85.94649168983183,84.92853925369353,84.86694214381284,84.62791851948175,84.48672796708584,84.16486234950084,83.94730150043614,83.88171428804769,83.82774277670238,83.55975235164334,83.41920860330092,82.59985442166085,82.36944496088982,82.05726016290478,77.3858228366621,76.57026961039973,76.12284670059285,75.84620494320025,75.32109008676477,75.04984429977227,74.81604018387954,74.58361820033605,74.56843162508834,74.51867261711911,74.35670041501308,74.19573987084722,74.14454483801123,73.44516234286183,73.12545115215933,73.11744465777583,73.05191117385225,71.97577532822777,71.49005471145058,71.37571241646155,71.28027873049473,70.80725840033709,70.8036729100543,70.59489107115465,70.5137912923328,70.32924850620888,70.10402260920615,69.87603030464756,69.53117718125154,68.76545788465243,68.7583591856931,68.11043541442514,67.84681011403126,67.57031378183422,67.49210331323822,66.93047730912309,66.71306083916912,66.39571825493401,66.22201681367976,66.19539050842751,65.99828309402207,65.90381340903588,65.81837592333682,65.72072397937451,65.45091124632094,65.2060580294703,65.18833104589127,64.91381911031287,64.01338306567423,63.55208562562408,63.408774691407906,62.96151654442945,62.15693255499719,61.704988716942594,61.68103571158919,61.12846369624528,60.94040505871836,60.53446681780395,59.60258979270896,59.501725048668156,58.742413649236255,57.49394383110114,54.32731105384171,54.29661161398992,53.970529417556335,53.469482771859724,52.70474076045028,51.38101859508939,50.140369305066365,49.40084938640584,48.442452387338925,48.38536052890897,45.84403019292141,45.47670698823559,45.41866897559957,45.41839217187456,45.29854611720612,45.211042509769136,45.15247606786069,44.938658094208186,44.77623598436978,44.729707097891634,44.62659634607659,42.56847210635154,42.50001401067301,42.38038561219,42.25147248296379,42.245345030701486,42.10547886717194,42.02962522632461,41.292287162752835,40.48369623343603,40.05187809911719,39.44116550297453,39.22623428954331,38.994746522073655,37.22262153464601,34.051162682640495,33.65516414161626,33.24546949495093,32.04444839294616,31.15093598738311,29.798721826132354,27.885322477001658,27.010938977679114,26.89767806161609,25.570356360313852,25.490008250305177,25.38273514075396,24.67652712609007,23.83479381888901,23.577190034623282,23.026230758386063,22.89963958324237,21.725530299425216,21.556648568603347,21.36027472173769,21.09349463960433,20.957866383366827,20.852247937344746,20.774202827642693,20.052440768387406,19.35892677385418,18.6977980136776,18.35686013324457,18.092108649551577,17.29334843585728,17.27460899582585,16.119283038362738,15.461985061347962,15.456634123225468,14.738922931615443,14.619956066410086,14.213845410890979,13.357781713932557,12.936051624355805,12.693293877613826,12.295631198551366,12.085153010540349,11.5622321506254,11.252575124466176,11.241104728298867,6.566564032023734,6.540356639138789,6.41422081661375,6.323430768076129,6.2725803662237425,6.18927334696072,6.087052810453401,5.988978382332227,5.953663700841609,5.948105629162042,4.899697738441464,3.344716682179103,2.8180585211394806,2.2537574560875893,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
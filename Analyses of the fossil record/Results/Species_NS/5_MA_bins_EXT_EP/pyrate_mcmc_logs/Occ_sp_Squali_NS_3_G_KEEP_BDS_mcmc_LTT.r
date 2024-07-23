

pdf(file='../Results/Species_NS/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_sp_Squali_NS_3_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_3_G_KEEP_BDS_mcmc"
ts=c(46.884367496229395, 63.516385646138644,67.13582605910057,84.39063049916918,85.90252831518816,18.86366590792056,28.9118255858142,6.944197835364086,83.95634511574465,78.87049579671636,6.487513206258119,68.99568416005854,84.16517496115891,22.53101060396193,22.451434048111217,27.25169578479279,53.35417790895794,58.2177256868321,69.37976430056176,84.34113979054285,48.78666174665638,30.777435495630566,69.91393015090325,20.740857452687255,58.72752010301787,66.93598565580491,33.51779964163651,6.711059617579662,23.259927682341974,49.428553423145935,78.24793666124283,69.80201101006124,76.62444440156962,101.50798369286608,106.65637515477377,109.278887571175,112.3124005149014,91.60131752462753,5.8331383968871275,67.1241898024014,49.05185222029314,46.77742909596384,34.06412572423655,76.67906178901934,6.406306450960904,66.98035112778864,33.27596127417251,69.36565630826186,68.82782650185025,65.53029119088033,12.363554194520724,2.381756270336479,82.2159823830262,53.54831072082508,53.34951636345449,77.96795042169266,21.846326815513237)
te=c(38.9886641183903, 61.332281593760555,63.55280589523365,65.26795846330424,83.90665834999966,0.0,17.33671834247834,0.0,61.34176909406666,64.93346472634323,0.0,65.40679434821625,71.43041501085372,0.0,0.0,19.25144069940385,46.49786442170514,48.58541758662022,65.38499844420394,69.12112256985341,40.315806151612755,14.03526790200108,62.16184279020072,3.6229176000302776,34.19810939899398,46.06579465899072,30.63277269220048,0.0,12.086006614606319,39.11519310807494,67.17667619509913,65.65643108961869,71.60999395219353,93.51339329823728,93.33325352040754,106.48378732709374,94.59451711387597,73.50836539368512,0.0,53.42460778453453,38.65095784355243,39.128741899096894,13.22746082569223,64.7321522922805,0.0,47.57873987896433,14.047105022074975,61.392634596968364,65.28779071330536,41.092254085488726,2.7121604528867778,0.0,70.22825918502954,39.337551646953614,39.277433530685286,72.20130192053925,0.0)
div_traj=c(0, 1,2,3,2,3,2,1,0,1,2,3,4,5,6,5,6,7,8,9,10,11,10,9,8,7,6,7,8,9,10,9,10,11,10,11,12,13,14,13,14,13,12,11,10,9,8,7,8,7,6,5,4,5,6,7,6,7,8,9,10,11,10,9,10,11,10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,3,4,5,6,7,8,9,10,9,10,9,8,7,6,7,6,7,8,9,10,11,10,9,10)
time_events=c(112.3124005149014, 109.278887571175,106.65637515477377,106.48378732709374,101.50798369286608,94.59451711387597,93.51339329823728,93.33325352040754,91.60131752462753,85.90252831518816,84.39063049916918,84.34113979054285,84.16517496115891,83.95634511574465,83.90665834999966,82.2159823830262,78.87049579671636,78.24793666124283,77.96795042169266,76.67906178901934,76.62444440156962,73.50836539368512,72.20130192053925,71.60999395219353,71.43041501085372,70.22825918502954,69.91393015090325,69.80201101006124,69.37976430056176,69.36565630826186,69.12112256985341,68.99568416005854,68.82782650185025,67.17667619509913,67.13582605910057,67.1241898024014,66.98035112778864,66.93598565580491,65.65643108961869,65.53029119088033,65.40679434821625,65.38499844420394,65.28779071330536,65.26795846330424,64.93346472634323,64.7321522922805,63.55280589523365,63.516385646138644,62.16184279020072,61.392634596968364,61.34176909406666,61.332281593760555,58.72752010301787,58.2177256868321,53.54831072082508,53.42460778453453,53.35417790895794,53.34951636345449,49.428553423145935,49.05185222029314,48.78666174665638,48.58541758662022,47.57873987896433,46.884367496229395,46.77742909596384,46.49786442170514,46.06579465899072,41.092254085488726,40.315806151612755,39.337551646953614,39.277433530685286,39.128741899096894,39.11519310807494,38.9886641183903,38.65095784355243,34.19810939899398,34.06412572423655,33.51779964163651,33.27596127417251,30.777435495630566,30.63277269220048,28.9118255858142,27.25169578479279,23.259927682341974,22.53101060396193,22.451434048111217,21.846326815513237,20.740857452687255,19.25144069940385,18.86366590792056,17.33671834247834,14.047105022074975,14.03526790200108,13.22746082569223,12.363554194520724,12.086006614606319,6.944197835364086,6.711059617579662,6.487513206258119,6.406306450960904,5.8331383968871275,3.6229176000302776,2.7121604528867778,2.381756270336479,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
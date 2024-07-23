

pdf(file='../Results/Genus/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_7_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_7_G_KEEP_BDS_mcmc"
ts=c(46.538053854880715, 48.66093605688474,65.55596553091759,85.72280791801924,68.12211432086993,86.08684873996894,38.87846421057016,88.2002906896011,85.56573650189794,61.285875992182184,63.75037384792548,31.094757045987446,64.57304204146854,85.43485698602574,49.828892814387835,46.85289993279957,47.142426291111036,25.57538559294674,75.18532383447045,86.56321604065654,57.06627693418615,68.28796560831515,67.93950833329144,50.735171664236184,36.055083888226854,48.9059855161477,75.40367259424195,97.66989732352357,128.20011140855297,94.27543936344784,81.35393745200726,49.08132058010719,35.338682105644125,40.093799045936194,97.88912366756178,48.043669470810265,100.71020751004531,48.68452942221206,34.935429828908646)
te=c(40.53429361411311, 11.164460412917958,60.9712462118806,66.34995306061333,0.0,0.0,0.0,68.20222540946912,65.34133831724839,0.0,0.0,19.918091272105293,47.00188728391576,65.5576910364837,14.730790738143872,0.0,0.0,0.0,70.15975447574364,61.18198096929758,0.0,51.24118736831382,66.2203268532726,11.963307097651608,0.0,10.324469491315458,66.84205167098975,67.65119363815812,67.42819992501576,66.71505458083978,0.0,0.0,0.0,0.0,12.640525457687856,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,14,13,14,15,14,13,12,11,10,9,8,9,8,9,10,11,10,9,10,9,10,11,12,13,14,15,16,17,16,17,18,17,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(128.20011140855297, 100.71020751004531,97.88912366756178,97.66989732352357,94.27543936344784,88.2002906896011,86.56321604065654,86.08684873996894,85.72280791801924,85.56573650189794,85.43485698602574,81.35393745200726,75.40367259424195,75.18532383447045,70.15975447574364,68.28796560831515,68.20222540946912,68.12211432086993,67.93950833329144,67.65119363815812,67.42819992501576,66.84205167098975,66.71505458083978,66.34995306061333,66.2203268532726,65.5576910364837,65.55596553091759,65.34133831724839,64.57304204146854,63.75037384792548,61.285875992182184,61.18198096929758,60.9712462118806,57.06627693418615,51.24118736831382,50.735171664236184,49.828892814387835,49.08132058010719,48.9059855161477,48.68452942221206,48.66093605688474,48.043669470810265,47.142426291111036,47.00188728391576,46.85289993279957,46.538053854880715,40.53429361411311,40.093799045936194,38.87846421057016,36.055083888226854,35.338682105644125,34.935429828908646,31.094757045987446,25.57538559294674,19.918091272105293,14.730790738143872,12.640525457687856,11.963307097651608,11.164460412917958,10.324469491315458,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
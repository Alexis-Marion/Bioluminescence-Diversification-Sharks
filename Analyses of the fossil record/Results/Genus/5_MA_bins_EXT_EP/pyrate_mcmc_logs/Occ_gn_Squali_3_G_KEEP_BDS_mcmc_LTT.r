

pdf(file='../Results/Genus/5_MA_bins_EXT_EP/pyrate_mcmc_logs/Occ_gn_Squali_3_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_3_G_KEEP_BDS_mcmc"
ts=c(46.699717984063476, 48.65847407231578,67.07729680414178,86.79629180927338,69.9586758680012,86.00085519782988,38.715968740943396,88.91312602011644,86.76386017403821,60.723323591029825,67.43077882052245,34.97961812974382,65.286327620935,84.94465907189966,48.305477690192454,46.87317176627648,47.18358753193905,24.663074737629497,78.34531706437346,87.32421040138975,60.10279649988852,68.78546799168754,68.01721382960535,57.94677611699125,34.328095538825515,48.71776773245422,76.96656278768916,96.5681432116688,128.99290105616316,93.58534548816664,81.56077345418905,48.78795868420895,33.08800963636994,41.570022385694216,96.27512043834679,46.846203965226735,97.47675400564638,49.03089051635223,35.84965515323878)
te=c(40.176607506140556, 12.062738478416602,61.59158110723355,65.22742111802165,0.0,0.0,0.0,70.13709803624582,65.52056492367491,0.0,0.0,17.497781977996024,46.25592796090049,64.92417097559165,14.746934226908687,0.0,0.0,0.0,73.91476149833352,62.6630638691445,0.0,47.752288971150925,65.5173817186289,10.564845027573991,0.0,11.614128384877054,64.92190205475505,69.5331867773367,65.62183943097291,65.3353569452929,0.0,0.0,0.0,0.0,15.591459543649327,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,12,13,12,13,14,15,16,15,14,13,12,13,12,11,10,9,8,9,10,11,12,13,14,15,16,15,16,17,18,19,18,19,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(128.99290105616316, 97.47675400564638,96.5681432116688,96.27512043834679,93.58534548816664,88.91312602011644,87.32421040138975,86.79629180927338,86.76386017403821,86.00085519782988,84.94465907189966,81.56077345418905,78.34531706437346,76.96656278768916,73.91476149833352,70.13709803624582,69.9586758680012,69.5331867773367,68.78546799168754,68.01721382960535,67.43077882052245,67.07729680414178,65.62183943097291,65.52056492367491,65.5173817186289,65.3353569452929,65.286327620935,65.22742111802165,64.92417097559165,64.92190205475505,62.6630638691445,61.59158110723355,60.723323591029825,60.10279649988852,57.94677611699125,49.03089051635223,48.78795868420895,48.71776773245422,48.65847407231578,48.305477690192454,47.752288971150925,47.18358753193905,46.87317176627648,46.846203965226735,46.699717984063476,46.25592796090049,41.570022385694216,40.176607506140556,38.715968740943396,35.84965515323878,34.97961812974382,34.328095538825515,33.08800963636994,24.663074737629497,17.497781977996024,15.591459543649327,14.746934226908687,12.062738478416602,11.614128384877054,10.564845027573991,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
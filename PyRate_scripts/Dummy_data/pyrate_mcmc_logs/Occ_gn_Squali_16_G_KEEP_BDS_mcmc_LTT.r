

pdf(file='../Results/Genus/5_MA_bins_EP/pyrate_mcmc_logs/Occ_gn_Squali_16_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_16_G_KEEP_BDS_mcmc"
ts=c(46.54063078692182, 48.60289653530871,66.35000777134185,87.50320084490924,70.22325934882562,86.98522099460793,39.34747798094595,89.43947692714642,85.4291505558703,58.8115068063998,66.12877855116433,33.59512908284807,65.4336146349955,85.51577392120592,48.32425390244977,46.81230610203543,47.15421774794787,20.600361067480492,77.44904065724862,86.84826123766312,59.92602264695125,68.11898879170299,73.37848526069365,51.2615656544196,38.09920071761691,49.11344980583705,77.08115664681816,97.00292546570496,130.99745312518712,94.12769592919085,83.81919250455543,48.523786207653835,35.20509909095625,41.41586483292331,96.71348399791043,46.74394023418966,100.12832518871635,49.062775625533746,33.7436878457173)
te=c(40.64434734589142, 10.214840584269997,60.059424394081155,65.66057161399542,0.0,0.0,0.0,68.72543310669374,66.78915962869526,0.0,0.0,17.060709776648817,47.18318848040716,66.72149948636383,15.08026685048832,0.0,0.0,0.0,72.31548450222803,57.988100138475254,0.0,50.07380431596107,68.78461355660613,11.451654634953648,0.0,11.584556743096156,65.75698047448184,68.79587450754792,66.64244494255588,66.98751921801038,0.0,0.0,0.0,0.0,15.351547124377293,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,15,14,13,12,13,12,11,10,9,10,11,10,9,10,9,10,11,10,11,10,11,12,13,14,15,14,15,16,17,18,19,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(130.99745312518712, 100.12832518871635,97.00292546570496,96.71348399791043,94.12769592919085,89.43947692714642,87.50320084490924,86.98522099460793,86.84826123766312,85.51577392120592,85.4291505558703,83.81919250455543,77.44904065724862,77.08115664681816,73.37848526069365,72.31548450222803,70.22325934882562,68.79587450754792,68.78461355660613,68.72543310669374,68.11898879170299,66.98751921801038,66.78915962869526,66.72149948636383,66.64244494255588,66.35000777134185,66.12877855116433,65.75698047448184,65.66057161399542,65.4336146349955,60.059424394081155,59.92602264695125,58.8115068063998,57.988100138475254,51.2615656544196,50.07380431596107,49.11344980583705,49.062775625533746,48.60289653530871,48.523786207653835,48.32425390244977,47.18318848040716,47.15421774794787,46.81230610203543,46.74394023418966,46.54063078692182,41.41586483292331,40.64434734589142,39.34747798094595,38.09920071761691,35.20509909095625,33.7436878457173,33.59512908284807,20.600361067480492,17.060709776648817,15.351547124377293,15.08026685048832,11.584556743096156,11.451654634953648,10.214840584269997,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                


pdf(file='../Results/Genus/5_MA_bins_EXT_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_5_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_5_G_KEEP_BDS_mcmc"
ts=c(46.34905558491144, 48.8584940120742,63.97535242523312,90.46756756615784,72.94930886206907,86.32083969318897,39.58115486206337,88.76529383670936,86.14951753464317,59.24126365787408,65.11632009302008,32.88822310346817,65.09997742087471,85.84034028941805,48.839996400453785,46.978842684375024,46.96044009510334,21.88445298103742,75.74667611324563,88.10575847603586,59.90575614074879,70.92584445207642,72.14442365705801,53.50227981621157,35.36025035340117,48.71919818618042,77.41337301372668,96.63989133083847,130.53823870440075,94.32968483605055,83.88026626843592,48.55022700702621,33.00415886665299,40.50052161176127,96.42903627168533,47.4293637741394,98.20044943407602,49.101406863513475,36.486402266071615)
te=c(39.93804298467545, 10.916265562821858,60.822554803014896,66.23650910859755,0.0,0.0,0.0,68.59394782204275,65.51925705023861,0.0,0.0,19.787952197650412,45.076218132583115,67.02629582670716,12.649417966252074,0.0,0.0,0.0,70.96386205619011,61.2654688768664,0.0,49.54684806864142,68.72143566990368,11.530881399947152,0.0,11.404257581402707,65.18084021211033,65.07197233620894,68.52447055106835,65.135979574059,0.0,0.0,0.0,0.0,12.933381929552446,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,15,16,15,14,13,12,11,10,9,8,9,10,9,10,9,8,9,10,11,10,11,12,13,14,15,16,17,18,19,18,19,18,19,20,21,22,23,24,23,22,21,20,19,18)
time_events=c(130.53823870440075, 98.20044943407602,96.63989133083847,96.42903627168533,94.32968483605055,90.46756756615784,88.76529383670936,88.10575847603586,86.32083969318897,86.14951753464317,85.84034028941805,83.88026626843592,77.41337301372668,75.74667611324563,72.94930886206907,72.14442365705801,70.96386205619011,70.92584445207642,68.72143566990368,68.59394782204275,68.52447055106835,67.02629582670716,66.23650910859755,65.51925705023861,65.18084021211033,65.135979574059,65.11632009302008,65.09997742087471,65.07197233620894,63.97535242523312,61.2654688768664,60.822554803014896,59.90575614074879,59.24126365787408,53.50227981621157,49.54684806864142,49.101406863513475,48.8584940120742,48.839996400453785,48.71919818618042,48.55022700702621,47.4293637741394,46.978842684375024,46.96044009510334,46.34905558491144,45.076218132583115,40.50052161176127,39.93804298467545,39.58115486206337,36.486402266071615,35.36025035340117,33.00415886665299,32.88822310346817,21.88445298103742,19.787952197650412,12.933381929552446,12.649417966252074,11.530881399947152,11.404257581402707,10.916265562821858,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
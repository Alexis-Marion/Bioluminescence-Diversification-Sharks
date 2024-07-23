

pdf(file='../Results/Species_NS/5_MA_bins_EXT_EP_FT/pyrate_mcmc_logs/Occ_sp_Squali_NS_5_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_sp_Squali_NS_5_G_KEEP_BDS_mcmc"
ts=c(46.53024983742151, 65.55470031818342,64.96981947401719,83.56787567917993,85.22173018055237,20.64441466420456,26.08785046509974,6.815869701708548,83.74097908627117,75.84814967397631,6.225846988449661,71.13905944587333,84.18291484710105,20.93607750100256,20.55209748381604,27.585356443166837,54.21160264753608,57.669862100165915,71.42583103646734,84.16847821951086,46.64653823702834,25.52121832363229,68.1815468543749,21.36262182011542,59.88005147299728,67.06169834575425,31.969809308462295,6.599282611299942,22.499360470941152,48.86135316254012,76.17319359091067,68.96523329944351,81.05846191250436,100.21986693393224,106.35041843858855,110.7304123255391,112.99845395284942,91.0701637075003,6.109719643777041,64.44755549588884,48.66031561564671,46.58686411412389,29.50314443112278,78.36216053667437,7.268641089832898,67.05921591638518,31.56382015261889,69.14694750497556,73.10663601692455,66.22036136560202,21.489030961604733,1.3900669190911192,82.26860580146608,54.09403847462504,54.306615013980924,76.72326015027629,23.8378700274877)
te=c(39.14358488417823, 62.77412779718619,61.11708926262532,66.35612774964275,83.19182851747739,0.0,13.417818407016817,0.0,61.393716172998104,65.15373156485796,0.0,65.55384357646486,70.61620360918015,0.0,0.0,20.981626987586793,48.069423856526804,48.026610871082845,67.98878282664963,70.18092320802894,40.24131617154974,20.383822723857477,64.66115990314397,1.6985597729753041,36.12041125964943,46.00898028424239,30.90930100386001,0.0,11.833007117302254,39.280576790044904,65.76192454618544,66.22612513428467,71.10631610066868,94.69855273056508,94.63708922100685,100.22713255214313,94.93618155521708,72.7989581424417,0.0,52.516930441312226,39.25195062800882,39.2884770018285,13.260195588579812,69.84689947243113,0.0,51.059191659958415,13.918061627324658,61.07698004378176,68.61911737875212,39.39695471579361,3.631443005057665,0.0,69.1530294470235,39.366295880523836,39.40109909162661,70.62009187365753,0.0)
div_traj=c(0, 1,2,3,2,3,2,1,0,1,2,3,4,5,6,5,6,7,8,9,10,11,12,11,12,13,12,11,10,9,8,7,8,9,8,9,8,9,10,9,8,9,8,9,8,7,8,7,8,7,6,5,4,5,6,7,8,9,8,7,8,9,8,7,8,9,10,9,8,7,6,5,4,3,2,1,0,1,2,1,2,3,4,5,6,7,8,9,8,9,10,11,10,9,8,7,6,7,8,9,10,11,10,9,10)
time_events=c(112.99845395284942, 110.7304123255391,106.35041843858855,100.22713255214313,100.21986693393224,94.93618155521708,94.69855273056508,94.63708922100685,91.0701637075003,85.22173018055237,84.18291484710105,84.16847821951086,83.74097908627117,83.56787567917993,83.19182851747739,82.26860580146608,81.05846191250436,78.36216053667437,76.72326015027629,76.17319359091067,75.84814967397631,73.10663601692455,72.7989581424417,71.42583103646734,71.13905944587333,71.10631610066868,70.62009187365753,70.61620360918015,70.18092320802894,69.84689947243113,69.1530294470235,69.14694750497556,68.96523329944351,68.61911737875212,68.1815468543749,67.98878282664963,67.06169834575425,67.05921591638518,66.35612774964275,66.22612513428467,66.22036136560202,65.76192454618544,65.55470031818342,65.55384357646486,65.15373156485796,64.96981947401719,64.66115990314397,64.44755549588884,62.77412779718619,61.393716172998104,61.11708926262532,61.07698004378176,59.88005147299728,57.669862100165915,54.306615013980924,54.21160264753608,54.09403847462504,52.516930441312226,51.059191659958415,48.86135316254012,48.66031561564671,48.069423856526804,48.026610871082845,46.64653823702834,46.58686411412389,46.53024983742151,46.00898028424239,40.24131617154974,39.40109909162661,39.39695471579361,39.366295880523836,39.2884770018285,39.280576790044904,39.25195062800882,39.14358488417823,36.12041125964943,31.969809308462295,31.56382015261889,30.90930100386001,29.50314443112278,27.585356443166837,26.08785046509974,25.52121832363229,23.8378700274877,22.499360470941152,21.489030961604733,21.36262182011542,20.981626987586793,20.93607750100256,20.64441466420456,20.55209748381604,20.383822723857477,13.918061627324658,13.417818407016817,13.260195588579812,11.833007117302254,7.268641089832898,6.815869701708548,6.599282611299942,6.225846988449661,6.109719643777041,3.631443005057665,1.6985597729753041,1.3900669190911192,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
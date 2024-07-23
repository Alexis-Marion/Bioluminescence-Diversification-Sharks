

pdf(file='../Results/Genus_NS/5_MA_bins_EP_FT/pyrate_mcmc_logs/Occ_gn_Squali_NS_19_G_KEEP_BDS_mcmc_LTT.pdf',width=0.6*20, height=0.6*10)
title = "Occ_gn_Squali_NS_19_G_KEEP_BDS_mcmc"
ts=c(49.237641493650834, 65.90313767069163,86.97005242697828,67.63627587380893,86.3227169420481,36.6790407228225,88.41495780534173,86.58513104769489,58.530969194359315,67.37999009681266,33.75886796753766,67.09104474981551,85.74945643746167,49.4568556464862,47.41253626863276,88.13810778548164,59.692205200492104,70.90728981863433,56.85687840222827,36.471985159334466,49.03321238978382,75.71159697993399,97.41733423992882,131.66070325762274,94.11390695038958,82.18305321791559,49.595922932696425,34.895526134919535,43.70330097460335,96.6378650938989,47.60540655090275,99.0704621537339,49.488254915516784,34.76004976272383)
te=c(9.281499775500116, 60.59438543960867,65.73586491407326,0.0,0.0,0.0,70.02215534387017,66.07635222361067,0.0,0.0,19.51887368197039,44.902370160776655,67.9526534407937,15.272877245455678,0.0,58.05572106218707,0.0,47.101066077649264,9.580666358286779,0.0,11.591819712900044,65.77583320328021,69.36643578047493,68.87229044237613,66.9750996945702,0.0,0.0,0.0,0.0,15.319030669637872,0.0,0.0,0.0,0.0)
div_traj=c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,13,12,11,10,11,12,13,12,11,12,11,10,9,10,11,10,11,12,13,14,15,16,17,18,17,16,17,18,19,20,21,22,21,20,19,18,17,16)
time_events=c(131.66070325762274, 99.0704621537339,97.41733423992882,96.6378650938989,94.11390695038958,88.41495780534173,88.13810778548164,86.97005242697828,86.58513104769489,86.3227169420481,85.74945643746167,82.18305321791559,75.71159697993399,70.90728981863433,70.02215534387017,69.36643578047493,68.87229044237613,67.9526534407937,67.63627587380893,67.37999009681266,67.09104474981551,66.9750996945702,66.07635222361067,65.90313767069163,65.77583320328021,65.73586491407326,60.59438543960867,59.692205200492104,58.530969194359315,58.05572106218707,56.85687840222827,49.595922932696425,49.488254915516784,49.4568556464862,49.237641493650834,49.03321238978382,47.60540655090275,47.41253626863276,47.101066077649264,44.902370160776655,43.70330097460335,36.6790407228225,36.471985159334466,34.895526134919535,34.76004976272383,33.75886796753766,19.51887368197039,15.319030669637872,15.272877245455678,11.591819712900044,9.580666358286779,9.281499775500116,0.0)
                par(mfrow=c(1,2))
                L = length(ts)
                plot(ts, 1:L , xlim=c(-max(ts)-1,0), pch=20, type="n", main=title,xlab="Time (Ma)",ylab="Lineages")
                for (i in 1:L){segments(x0=-te[i],y0=i,x1=-ts[i],y1=i)}    
                t = -time_events
                plot(div_traj ~ t,type="s", main = "Diversity trajectory",xlab="Time (Ma)",ylab="Number of lineages",xlim=c(-max(ts)-1,0))
                abline(v=c(65,200,251,367,445),lty=2,col="gray")
                
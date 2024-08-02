 # Ecological opportunity fostered deep-sea  diversification of squaliform sharks

## Summary 

- [Summary](#Summary)
- [Overview](#Overview)
- [1 Bayesian estimation of deep-time diversfiication with PyRate](#1-Bayesian-estimation-of-deep-time-diversification-with-PyRate)
	- [1.1 Preservation model](#11-Preservation-model)
	- [1.2 Running PyRate](#12-Running-Pyrate)
	- [1.3 Assess convergence](#13-Assess-convergence)
  	- [1.4 Plotting PyRate results](#14-Plotting-PyRate-results)
- [2 Analyses of the fossil record](#2-Analyses-of-the-fossil-record)
  - [2.1 Selecting PyRate model](#21-Selecting-PyRate-model)
	- [2.2 Extracting time for speciation and extinction](#22-Extracting-time-for-speciation-and-extinction)
	- [2.3 Estimating lineage through time per habitat](#23-Estimating-lineage-through-time-per-habitat)
  	- [2.4 Estimating tempo of origination](#24-Estimating-tempo-of-origination)
	- [2.5 Grafting fossils](#25-Grafting-fossils)
- [3 Phylogenetic comparative analyses](#3-Phylogenetic-comparative-analyses)
	- [3.1 Analyses of discrete trait evolution with corHMM](#31-Analyses-of-discrete-trait-evoplution-with-corHMM)
 	- [3.2 Analyses of continuous trait evolution with OUwie and phylogenetic ANOVA](#32-Analyses-of-continuous-trait-evolution-with-OUwie-and-phylogenetic-ANOVA)
  	- [3.3 Joint estimation of discrete and continuous trait evolution with hOUwie](#32-Joint-estimation-of-discrete-and-continuous-trait-evolution-with-hOUwie)
- [Reference](#Reference)

<p align="justify"> This repository's purpose is to give a means of replicability to the article "Ecological opportunity fostered deep-sea  diversification of squaliform sharks" but can be generalized to other similar data. All of the presented scripts are written in R language (R Core Team, 2022). Each script is available as both an annotated notebook (.ipynb) or a raw .r file (unannotated).
	If you plan to use any of these scripts, please cite "XXX". </p>

## Overview

<p align="justify"> This repository contains scripts and html files for performing the following analyses:

**1**: Bayesian estimation of deep-time diversification with PyRate

**2**: Analyses of the fossil record

**3**: Phylogenetic comparative analyses

<p align="justify"> All data used to perform each analyses are deposited on this repository </p>

## 1 Bayesian estimation of deep-time diversification with PyRate

`package requirement (corHMM, mclust, string, phytools, qpcR, ggtree, secsse, ggtree, phytools, treedataverse, RColorBrewer, ggplot2, ggpmisc, optional(stringr))`

`used script (corHMM.r; corHMM_Diet.r; Plot_ASE.r)`

<p align="justify">  In this first session, we will be using PyRate (Silvestro et al, 2014). PyRate is a programm implemented in Python whose aiml is to jointly estimate the preservation process, the tempo of origination and extinction of lineages based on their occurences in the fossil record. Here, we will assume that the PyRate repository with its functions is at the root of the current working directory.</p>

### 1.1 Preservation model

<p align="justify"> One the main strength of PyRate is its ability to account for the biais of the fossil record, by estimating a preservation process correcting the estimated age derived from raw occurences data. Thus, choosing the best-fit preservation model for any PyRate analysis is critical. Fortunately, Silvestro et al. (2018) implemented a likelihood-based approach for preservation model selection.Yet, while this procedure is certainly useful, it is incomplete. Indeed the first implementation allowed for model selection across HPP, NHPP, TPP and alternative version of the TPP, with missing bins. Howver, bin removal occured only once, and were not recursive. Consequently, model selection is incomplete. Here, we corrected and enhanced this procedure, by performing model selection on all PyRate replicate (here 100). Furthermore, we allowed for recursivec bin removal, meaning that the best fit TPP model could be a two-bin model whereas the generating TPP model could be a five-bin model. Model selection is performed with pairwise comparaison of the AICc metrics across all replicates. </p>

</p> 

### 1.2 Running PyRate

<p align="justify"> The script provided in this section is rather simple, and run a BDCS analyses on 20,000,000 generation on the genus dataset including singeltons, with diversification shifts every 5 Myrs and integrating preservation shifts from the 1.1 section. Here, to ber computaionally efecient, we choose to parallelise our run on 20 cpu.</p> 

### 1.3 Assess convergence

PyRate is a bayesian programm, thus we will consider that a PyRate run is finished when it achieved convergence. A popular metric to evaluate convergence is the ESS (effective sample size), and it is generally considered that when its number is above 200, convergence is achevied. Thus, we assessed convergence on all run using the "" scripts. Furthermore, these scritps give additional usefull metrics on the run, such as the origination age, extinction age (if including solely extinct taxa),and the proportion of Ts and Te, with ESS above 200. We also provided a graphical output that can be executed using ....

### 1.4 Plotting PyRate results

In this last section, we provided plotting scripts to display graphically each PyRate output. These scripts will take as input the output direcotry of a regular (BDCS or RJMCMC) PyRate run, and will represent, the RTT (origination and extiction), the diversification RTT, the LTT and the QTT.

## 2 Diversification analyses

`package requirement (ape, mclust, secsse, DDD, tidyverse, qgraph, ggtree, phytools, treedataverse, RColorBrewer, ggplot2, ggpmisc, optional(stringr))`

`used script (SecSSE_Size.r; SecSSE_Reproduction.r; SecSSE_Habitat.r; SecSSE_Diet.r; SecSSE_ASE.r; Plot_ASE.r)`

<p align="justify"> Now that we pictured trait evolution dynamics, we may want to assess whether our examined traits are responsible for extant diversity patterns. To do so we rely on SSE models (State-dependent speciation and extinction models), which are well-known for accounting for the impact that trait evolution has on patterns of lineage diversification. However, accounting for trait-dependent diversification is subject to numerous methodological biases (Beaulieu and Donoghue, 2013). Indeed, SSE models can falsely indicate an effect of the focal trait on diversification. Models with hidden traits (aka concealed traits), such as SeCSSE (Herrera-Alsina et al., 2019) or HiSSE (Beaulieu and O'Meara, 2016) can account for hidden variables in trait-dependant diversification. Thus, we will be using SecSSE, an SSE implementation for detecting trait-dependent diversification using phylogenies on multi-state characters, while being robust to false-postive. </p>

### 2.1 Data cleaning and model construction

<p align="justify"> Highly similar to corHMM, we will first need to discretise our continuous traits and clean our data. However, unlike corHMM, SecSSE estimates jointly transition rates and diversification rates (speciation and extinction). This is both a good thing, as it has been shown that accounting for diversification allows for better estimation of trait evolution (Maddison, 2006), and a bad thing as it dramatically increases the number of parameters in our analyses. Over-parametrization is a common issue in statistics, but it has not often been addressed in macroevolution studies. Here, we want to avoid this as much as possible by limiting the upper number of estimated parameters. To do so, we kept the structure of the best-fitting corHMM model for each trait, while still estimating its transition rates. Namely, if the best-fitting model for reproduction is a two-rate sequential symmetrical model, we will forbid direct transition from A to C, while constraining transitions from A to B to be equal to B to A. To minimize the effect of structuration on the output of the model, we kept the same structure of the transition matrix across all variants for each trait.

For each trait, we constructed seven models, one for constant rates (CR), three for examined-trait diversification (ETD) and three for concealed-trait diversification (CTD). The three variants for CTD and ETD models include a pure speciation model (different speciation for each state, one shared extinction), a pure extinction model (different extinction for each state, one shared speciation) and a net diversification model (different speciation and extinction for each state). For each of the seven models to avoid finding a local optimum, following Herrera-Alsina et al. (2019), we used three sets of initial parameters. One set was estimated using the standard birth-death model (bd_ML function in the R package “DDD” 5.2.2; Etienne et al., 2012), one was halved, and the other was doubled.
</p>

### 2.2 Running and selecting the likeliest model

<p align="justify"> Then, we may perform each model and save relevant metrics for comparison. Here, we save the output of each of the 21 models for each trait in a separate directory. SecSSE does not compute the AICc directly, thus we used the following formula to estimate the AICc "2*(-(ll)+2*k+(2*k*(k+1))/(n-k-1))" where ll is the log-likelihood of the model, k its number of parameters and n the number of taxa included in the analysis. Then directly in the script, we construct a data frame filled with the log-likelihood, the number of parameters, the AICc, the $\Delta$ AICc, the $\omega$ AICc and the estimated rates. </p>


### 2.3 Ancestral state estimation

<p align="justify"> After the selection of the best-fitting model, one may want to estimate ancestral states across the phylogeny. However, unlike corHMM, SecSSE does not directly infer ancestral state, rather we may require another procedure to properly estimate ancestral state. This procedure is highlighted in the script "SecSSE_ASE.r", where there is the code for filling proper parameters. This procedure is a bit tricky, but you will find in the code that everything is annotated. Finally, after running our model (which is quite fast), we may plot our ancestral state directly on the phylogeny, similarly to what was done with corHMM  </p>

## 3 Sensitivity analyses 

`package requirement (ape, mclust, secsse, DDD, tidyverse, qgraph, tidyverse, ggpubr, rstatix, optional(string))`

`used script (corHMM.r; corHMM_Diet.r; SecSSE_Size.r; SecSSE_Reproduction.r; SecSSE_Habitat.r; SecSSE_Diet.r; Posterior_test_comparative_analysis.r)`

<p align="justify"> Both trait-dependent evolution and diversification model require phylogeny. However, phylogenetic trees are hypotheses. Consequently, a phylogeny may not accurately reflect, which may affect downstream analyses. Worse, when a calibrated phylogeny is employed the uncertainty is twofold: as both dating and phylogenetic estimation are affected. We accounted for both phylogenetic and dating uncertainty by performing sensitivity analyses on 100 trees extracted from the posterior distribution of the BEAST analysis. To perform these analyses, the user just has to replace the input tree and the data frame output manually in each script (for an example of how to do it, see the following scripts "corHMM_replicated.r", "SecSSE_Size_replicated.r", and "multiRun_script.sh"). We analyzed whether we could reliably recover the best-fitting model for each trait and analysis (corHMM, SecSSE) by performing a non-parametric alternative of the repeated measure ANOVA (Friedman test) comparing the AICc of each model and then performing pairwise comparisons of each model with a paired signed-rank Wilcoxon test using the R package “rstatix” (Kassambara, 2023) where we reported the T statistic and p-value. These tests are implemented in the script "Posterior_test_comparative_analysis.r", which will require a directory with all data frame replicates (either corHMM or SecSSE) as input.</p>


### Reference

Beaulieu, J. M. Donoghue, M. J. (2013). Fruit evolution and diversification in campanulid angiosperms: Campanulid fruit evolution. Evolution. 67(11): 3132-3144.

Beaulieu, J. M., & O’Meara, B. C. (2016). Detecting Hidden Diversification Shifts in Models of Trait-Dependent Speciation and Extinction. Systematic Biology, 65(4), 583-601. 

Boyko, J. D., & Beaulieu, J. M. (2021). Generalized hidden Markov models for phylogenetic comparative datasets. Methods in Ecology and Evolution, 12(3), 468-478. 

Etienne, R. S., Haegeman, B., Stadler, T., Aze, T., Pearson, P. N., Purvis, A., & Phillimore, A. B. (2012). Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record. Proceedings of the Royal Society B: Biological Sciences, 279(1732), 1300-1309. 

Herrera-Alsina, L. van Els, P. Etienne, R. S. (2019). Detecting the dependence of diversification on multiple traits from phylogenetic trees and trait data. Systematic Biology. 68(2): 317-328.

Fraley, C. Raftery, A. E. Murphy, T. B. & Scrucca, L. (2012). mclust Version 4 for R: normal mixture modelling for model-based clustering, classification, and density estimation. Technical Report No. 597. Seattle, WA: Department of Statistics, University of Washington.

Maddison, W. P. (2006). Confounding asymmetries in evolutionary diversification and character change. Evolution, 60(8), 1743-1746.

R Core Team (2022). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

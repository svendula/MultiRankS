# Estimating the underlying signal from multiple ranked lists
Vendula Svendova, Michael G. Schimek

Application example from:

Fortney et al.(2015). Prioritizing therapeutics for lung cancer: an integrative meta-analysis of cancer gene signatures and chemogenomic data. PLoS Computational Biology.11(3).

We used a subset of 50 drugs, out of which 10 were ranked at the top of the 247 significant drugs chosen by Fortney et al (2015) as those consistently reversing lung cancer gene changes. The rest 40 drugs were chosen randomly.

Licence terms: http://creativecommons.org/licenses/by/4.0/

Fortney50_input.csv - Input rank matrix 

Fortney50_single_est.csv - The signal estimate from the Input rank matrix

Fortney50_boot.mean_est.csv - Mean of 30 bootstrap estimates (robust to contamination by 'bad exoerts')

Fortney50_groundtruth_median.csv - Median of the Connectivity Scores (ground truth)

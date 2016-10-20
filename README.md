# Estimating the underlying signal from multiple ranked lists
Vendula Svendova, Michael G. Schimek

Here the source code and genomic application data and results are provided.

## Source code

[`MultiRankS_funs.R`](https://github.com/svendula/Estimating-the-underlying-signal-from-multiple-ranked-lists/blob/master/MultiRankS_funs.R): The source code of the algorithm.

[`README_MultiRankS`](https://github.com/svendula/Estimating-the-underlying-signal-from-multiple-ranked-lists/blob/master/README_MultiRankS.md): a description of all the functions in [`MultiRankS_funs.R`](https://github.com/svendula/Estimating-the-underlying-signal-from-multiple-ranked-lists/blob/master/MultiRankS_funs.R)


## Application example

The application example was based on data from:

Fortney et al. (2015). _Prioritizing therapeutics for lung cancer: an integrative meta-analysis of cancer gene signatures and chemogenomic data_. PLoS Computational Biology. 11(3). <http://dx.doi.org/10.1371/journal.pcbi.1004068>

We used a subset of 30 drugs, out of which 10 were ranked at the top of the 247 most significant drugs chosen by Fortney et al. (2015) as those consistently reversing lung cancer gene changes. The remaining 20 drugs were chosen randomly.

The data from Fortney et al. were provided under the terms of the Creative Commons Attribution (CC BY) license, allowing for the data's reuse, redistribution, and adaption, providing the authors and original source are cited. For full terms see <http://creativecommons.org/licenses/by/4.0/>.

Fortney, K., Griesman, J., Kotlyar, M., Pastrello, C., Angeli, M., Sound-Tsao, M., and Jurisica, I. (2015). Prioritizing therapeutics for lung cancer: an integrative meta-analysis of cancer gene signatures and chemogenomic data. PLoS Computational Biology, 11(3):e1004068.


[`Fortney30_input.csv`](https://github.com/svendula/Estimating-the-underlying-signal-from-multiple-ranked-lists/blob/master/Fortney50_input.csv): Input rank matrix of the application example

[`Fortney30_groundtruth_median.csv`](https://github.com/svendula/Estimating-the-underlying-signal-from-multiple-ranked-lists/blob/master/Fortney50_groundtruth_median.csv): Median of the Connectivity Scores (ground truth) from the application example







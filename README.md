# MultiRankS

MultiRankS is an R tool for estimating the underlying signal from multiple ranked lists. It takes a matrix of ranks as an input, where rankers are in columns and ranked objects in rows. It outputs an estimated real signal value for each object. These signals are responsible for the rankings in a sense that they are the true values that the rankers intended to rank. The tool also calculates the standard errors of the signal estimates using either Metropolis-Hastings or adaptive MCMC optimisation.

For more details of the method, please refer to

Švendová, Vendula, and Michael G. Schimek (2017). A novel method for estimating the common signals for consensus across multiple ranked lists. *Computational Statistics & Data Analysis*, Volume 115, pp 122-135.
http://www.sciencedirect.com/science/article/pii/S0167947317301056

Here the source code and genomic application data and results are provided, under the terms of the LGPLv3 [`license`](https://github.com/svendula/MultiRankS/blob/master/LICENSE).

## Source code
Currently, the code cannot be installed directly from GitHub. To use the method, you need to save the file `MultiRankS_funs.R` locally and run the script in R using `source('MultiRankS_funs.R')`.

[`MultiRankS_funs.R`](https://github.com/svendula/MultiRankS/blob/master/MultiRankS_funs.R): The source code of the algorithm.

[`README_MultiRankS`](https://github.com/svendula/MultiRankS/blob/master/README_MultiRankS.md): a description of all the functions in [`MultiRankS_funs.R`](https://github.com/svendula/MultiRankS/blob/master/MultiRankS_funs.R)


## Example of usage
Set up parameters (generally, larger values add to precision but also computatinal time)
```r
source('MultiRankS_funs.R')
p = 5 # number of objects
n = 10 # number of rankers
l_max = 2 # maximum height of the sliding window (larger value = longer computation)
num.boot = 10 # number of bootstrap samples 
chain.length = 1000 # length of the MCMC chain
num.sim = 5 # number of independent chains
```
Create a random matrix of ranks `R.input` - objects in rows, rankers in columns:
```r
R.input = generate.random.rank.matrix(p,n)
```
Create bootstrap matrices
```r
boots = generate.bootstrap.samples(R.input,num.boot)
```
Calculate the list of probability matrices for the input:
```r
F.input = F_perm(R.input, l_max) 
```
Run adaptive MCMC
```r
results.by.boot = run.adaptiveMCMC.example(F.input, boots, num.sim, chain.length)
estimate = gather.results(results.by.boot)
```

## Application example in the reference article

The application example in the paper above was based on data from:

Fortney, K., Griesman, J., Kotlyar, M., Pastrello, C., Angeli, M., Sound-Tsao, M., and Jurisica, I. (2015). Prioritizing therapeutics for lung cancer: an integrative meta-analysis of cancer gene signatures and chemogenomic data. *PLoS Computational Biology*, 11(3):e1004068. <http://dx.doi.org/10.1371/journal.pcbi.1004068>

We used a subset of 30 drugs, out of which 10 were ranked at the top of the 247 most significant drugs chosen by Fortney et al. (2015) as those consistently reversing lung cancer gene changes. The remaining 20 drugs were chosen randomly.

The data from Fortney et al. were provided under the terms of the Creative Commons Attribution (CC BY) license, allowing for the data's reuse, redistribution, and adaption, providing the authors and original source are cited. For full terms see <http://creativecommons.org/licenses/by/4.0/>.


[`Fortney30_input.csv`](https://github.com/svendula/MultiRankS/blob/master/Fortney30_input.csv): Input rank matrix of the application example

[`Fortney30_groundtruth_median.csv`](https://github.com/svendula/MultiRankS/blob/master/Fortney30_groundtruth_median.csv): Median of the Connectivity Scores (ground truth) from the application example







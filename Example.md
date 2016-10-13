MultiRankS example
================
Vendula Svendova
5 September 2016

This is an example of the MultiRankS algorithm, as described in paper *Svendova, Schimek: A novel method for estimating the underlying signal from multiple ranked lists*. Here the number of objects, as well as assessors, is 5.

Generate the input
------------------

The assumed model is \[X.\text{input}_{ij} = \text{theta.true}_{i} + Z_{ij},\] where \(i=1,\ldots,p, j=1,\ldots,n\), and \(Z\sim N(0,\sigma_j^2)\)

### Simulate the true underlying signal `theta.true`

In real-world applications, `theta.true` is not known. Here, we simulate it in order to evaluate our estimate. We aim to estimate the normalised true signal `theta.true.norm` and its ranking `true.rank`.

``` r
source("MultiRankS_funs.R")
p = 10  # number of objects 
n = 10  # number of assessors 

set.seed(123)
theta.true = sort(rnorm(p, 0, 1), decreasing = TRUE)  # sorted true signal 
theta.true.norm = theta.true/norm_vec(theta.true)  # normalised true signal
theta.true.norm
```

    ##  [1]  0.59736152  0.54290209  0.16053829  0.04503125  0.02455825
    ##  [6] -0.08017141 -0.15522520 -0.19521510 -0.23923261 -0.44062407

``` r
true.rank = rank(-theta.true.norm)  # the true rank - the largest value will have rank 1, the smallest rank 5
true.rank
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

``` r
plot(theta.true.norm, col = "red", ylim = c(-1, 1), type = "b", 
    pch = 16, xlab = "Object", ylab = "Normalised true signal")
```

![](Example_markdown_files/figure-markdown_github/unnamed-chunk-1-1.png)

### Simulate the measurements `X.input` .

In real-world examples, `X.input` is either not known or contains incomparable values (e.g. one column has real values between 0 and 1, another has integers between 1 and 100). In order to simulate `X.input`, one has to choose a reasonable variance added random variable \(Z\sim N(0,\sigma_j^2)\). In this example, we choose each \(\sigma_j\) as the absolute random value from normal distribution with 0 mean and standard deviation 0.3.

``` r
set.seed(123)
sigmas = abs(rnorm(n, 0, 0.4))  # random standard deviations for each assessors 
X.input = matrix(nrow = p, ncol = n)  # the matrix of observed measurements 
for (i in 1:n) {
    set.seed(i)
    X.input[, i] = theta.true + rnorm(p, 0, sigmas[i])
}
X.input  # objects in rows, assessors in columns
```

    ##             [,1]        [,2]       [,3]       [,4]       [,5]        [,6]
    ##  [1,]  1.5746201  1.63248517  1.1153155  1.7211782  1.6715801  1.90002170
    ##  [2,]  1.5998794  1.57572756  1.3763234  1.5434082  1.6303006  1.12652194
    ##  [3,]  0.2735764  0.60711071  0.6222663  0.4860495  0.3959883  1.05683943
    ##  [4,]  0.4869342  0.02521292 -0.5890473  0.1460964  0.1329152  1.31418876
    ##  [5,]  0.1443808  0.06311953  0.1925757  0.1166383  0.1590157  0.08710174
    ##  [6,] -0.4141185 -0.21798542 -0.2113957 -0.2107376 -0.2613569  0.02229735
    ##  [7,] -0.3363851 -0.38047987 -0.3924054 -0.4817974 -0.4700801 -1.34381015
    ##  [8,] -0.3949504 -0.58254488  0.1357122 -0.5664870 -0.5933339 -0.05376180
    ##  [9,] -0.5577683 -0.50414036 -1.4467901 -0.6333641 -0.7016317 -0.65606882
    ## [10,] -1.3335263 -1.27783949 -0.4748780 -1.2149477 -1.2579190 -1.98428897
    ##              [,7]       [,8]        [,9]       [,10]
    ##  [1,]  2.13675670  1.6722624  1.50439457  1.71840677
    ##  [2,]  1.33806373  1.9839714  1.33439362  1.52586257
    ##  [3,]  0.33291194  0.2263826  0.42203067  0.21645626
    ##  [4,]  0.05327473 -0.1494483  0.05301821  0.02247723
    ##  [5,] -0.10845124  0.4429629  0.19037985  0.12301542
    ##  [6,] -0.40482416 -0.2847681 -0.55626020 -0.16069089
    ##  [7,] -0.30773015 -0.5318325 -0.11817413 -0.66101941
    ##  [8,] -0.58203827 -1.1111982 -0.56547328 -0.62530627
    ##  [9,] -0.65870790 -2.2105188 -0.75501190 -0.97683131
    ## [10,] -0.86130267 -1.5652220 -1.36477493 -1.31078230

### Construct the input rank matrix `R.input`.

The rank matrix `R.input` is the only input for our algorithm. Its columns are ranked columns of `X.input`.

``` r
R.input = matrix(nrow = p, ncol = n)  # the observed rankings 
R.input = apply(X.input, 2, function(x) rank(-x))
R.input  # objects in rows, assessors in columns
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,]    2    1    2    1    1    1    1    2    1     1
    ##  [2,]    1    2    1    2    2    3    2    1    2     2
    ##  [3,]    4    3    3    3    3    4    3    4    3     3
    ##  [4,]    3    5    9    4    5    2    4    5    5     5
    ##  [5,]    5    4    4    5    4    5    5    3    4     4
    ##  [6,]    8    6    6    6    6    6    7    6    7     6
    ##  [7,]    6    7    7    7    7    9    6    7    6     8
    ##  [8,]    7    9    5    8    8    7    8    8    8     7
    ##  [9,]    9    8   10    9    9    8    9   10    9     9
    ## [10,]   10   10    8   10   10   10   10    9   10    10

``` r
cor(R.input, method = "spearman")  # correlation between the assessors
```

    ##            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
    ##  [1,] 1.0000000 0.8909091 0.6848485 0.9393939 0.9151515 0.8787879
    ##  [2,] 0.8909091 1.0000000 0.7454545 0.9757576 0.9878788 0.8787879
    ##  [3,] 0.6848485 0.7454545 1.0000000 0.7454545 0.8060606 0.5636364
    ##  [4,] 0.9393939 0.9757576 0.7454545 1.0000000 0.9878788 0.9272727
    ##  [5,] 0.9151515 0.9878788 0.8060606 0.9878788 1.0000000 0.8909091
    ##  [6,] 0.8787879 0.8787879 0.5636364 0.9272727 0.8909091 1.0000000
    ##  [7,] 0.9636364 0.9636364 0.7333333 0.9878788 0.9757576 0.8909091
    ##  [8,] 0.9030303 0.9393939 0.8303030 0.9393939 0.9636364 0.8303030
    ##  [9,] 0.9393939 0.9757576 0.7939394 0.9757576 0.9878788 0.8545455
    ## [10,] 0.9030303 0.9636364 0.8303030 0.9757576 0.9878788 0.9151515
    ##            [,7]      [,8]      [,9]     [,10]
    ##  [1,] 0.9636364 0.9030303 0.9393939 0.9030303
    ##  [2,] 0.9636364 0.9393939 0.9757576 0.9636364
    ##  [3,] 0.7333333 0.8303030 0.7939394 0.8303030
    ##  [4,] 0.9878788 0.9393939 0.9757576 0.9757576
    ##  [5,] 0.9757576 0.9636364 0.9878788 0.9878788
    ##  [6,] 0.8909091 0.8303030 0.8545455 0.9151515
    ##  [7,] 1.0000000 0.9272727 0.9878788 0.9515152
    ##  [8,] 0.9272727 1.0000000 0.9515152 0.9515152
    ##  [9,] 0.9878788 0.9515152 1.0000000 0.9636364
    ## [10,] 0.9515152 0.9515152 0.9636364 1.0000000

### Calculate the list of probability matrices `F.input`

For each window size \(\ell=1,\ldots,\ell_0\), a \((p-\ell+1)\times p^{\ell}\) probability matrix is calculated and saved to `F.input`.

``` r
l_0 = 2  # maximum window size 
F.input = F_perm(R.input, l_0)
F.input
```

    ## $`l=1`
    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,]  0.7  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0     1
    ##  [2,]  0.3  0.9  1.0  1.0  1.0  1.0  1.0  1.0  1.0     1
    ##  [3,]  0.0  0.0  0.7  1.0  1.0  1.0  1.0  1.0  1.0     1
    ##  [4,]  0.0  0.1  0.2  0.4  0.9  0.9  0.9  0.9  1.0     1
    ##  [5,]  0.0  0.0  0.1  0.6  1.0  1.0  1.0  1.0  1.0     1
    ##  [6,]  0.0  0.0  0.0  0.0  0.0  0.7  0.9  1.0  1.0     1
    ##  [7,]  0.0  0.0  0.0  0.0  0.0  0.3  0.8  0.9  1.0     1
    ##  [8,]  0.0  0.0  0.0  0.0  0.1  0.1  0.4  0.9  1.0     1
    ##  [9,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.2  0.8     1
    ## [10,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.1  0.2     1
    ## 
    ## $`l=2`
    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ##  [1,]    0  0.3  0.3  0.3  0.3  0.3  0.3  0.3  0.3   0.3   0.6   0.9   0.9
    ##  [2,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [3,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [4,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [5,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [6,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [7,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [8,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##  [9,]    0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0   0.0   0.0   0.0   0.0
    ##       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
    ##  [1,]   0.9   0.9   0.9   0.9   0.9   0.9   0.9   0.7   1.0   1.0   1.0
    ##  [2,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.1   0.7   0.7   0.7
    ##  [3,]   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.0   0.0   0.0   0.2
    ##  [4,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [5,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [6,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [7,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35]
    ##  [1,]   1.0   1.0   1.0   1.0   1.0   1.0   0.7   1.0   1.0   1.0   1.0
    ##  [2,]   0.7   0.7   0.7   0.7   0.7   0.7   0.3   0.9   1.0   1.0   1.0
    ##  [3,]   0.2   0.2   0.2   0.2   0.2   0.2   0.0   0.0   0.2   0.4   0.4
    ##  [4,]   0.1   0.1   0.1   0.1   0.1   0.1   0.0   0.0   0.0   0.0   0.5
    ##  [5,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [6,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [7,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46]
    ##  [1,]   1.0   1.0   1.0   1.0   1.0   0.7   1.0   1.0   1.0   1.0   1.0
    ##  [2,]   1.0   1.0   1.0   1.0   1.0   0.3   0.9   1.0   1.0   1.0   1.0
    ##  [3,]   0.4   0.4   0.4   0.4   0.4   0.0   0.0   0.6   0.9   0.9   0.9
    ##  [4,]   0.5   0.5   0.5   0.6   0.6   0.0   0.1   0.2   0.4   0.9   0.9
    ##  [5,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [6,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [7,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57]
    ##  [1,]   1.0   1.0   1.0   1.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [2,]   1.0   1.0   1.0   1.0   0.3   0.9   1.0   1.0   1.0   1.0   1.0
    ##  [3,]   0.9   0.9   0.9   0.9   0.0   0.0   0.6   0.9   0.9   0.9   0.9
    ##  [4,]   0.9   0.9   1.0   1.0   0.0   0.1   0.2   0.4   0.9   0.9   0.9
    ##  [5,]   0.0   0.0   0.0   0.0   0.0   0.0   0.1   0.5   0.7   0.7   0.7
    ##  [6,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.2
    ##  [7,]   0.1   0.1   0.1   0.1   0.0   0.0   0.0   0.0   0.0   0.0   0.1
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66] [,67] [,68]
    ##  [1,]   1.0   1.0   1.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [2,]   1.0   1.0   1.0   0.3   0.9   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [3,]   0.9   0.9   0.9   0.0   0.0   0.6   0.9   0.9   0.9   0.9   0.9
    ##  [4,]   0.9   1.0   1.0   0.0   0.1   0.2   0.4   0.9   0.9   0.9   0.9
    ##  [5,]   0.7   0.7   0.7   0.0   0.0   0.1   0.6   0.9   0.9   0.9   0.9
    ##  [6,]   0.3   0.3   0.3   0.0   0.0   0.0   0.0   0.0   0.5   0.7   0.8
    ##  [7,]   0.1   0.1   0.1   0.0   0.0   0.0   0.0   0.0   0.1   0.2   0.3
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,69] [,70] [,71] [,72] [,73] [,74] [,75] [,76] [,77] [,78] [,79]
    ##  [1,]   1.0   1.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [2,]   1.0   1.0   0.3   0.9   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [3,]   0.9   0.9   0.0   0.0   0.6   0.9   0.9   0.9   0.9   0.9   0.9
    ##  [4,]   1.0   1.0   0.0   0.1   0.2   0.4   0.9   0.9   0.9   0.9   1.0
    ##  [5,]   0.9   0.9   0.0   0.0   0.1   0.6   1.0   1.0   1.0   1.0   1.0
    ##  [6,]   0.8   0.8   0.0   0.0   0.0   0.0   0.0   0.6   0.8   0.9   0.9
    ##  [7,]   0.4   0.4   0.0   0.0   0.0   0.0   0.0   0.3   0.7   0.8   0.9
    ##  [8,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.1   0.1   0.2
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
    ##       [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90]
    ##  [1,]   1.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [2,]   1.0   0.3   0.9   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [3,]   0.9   0.0   0.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [4,]   1.0   0.0   0.1   0.2   0.4   0.9   0.9   0.9   0.9   1.0   1.0
    ##  [5,]   1.0   0.0   0.0   0.1   0.6   1.0   1.0   1.0   1.0   1.0   1.0
    ##  [6,]   0.9   0.0   0.0   0.0   0.0   0.0   0.7   0.9   1.0   1.0   1.0
    ##  [7,]   0.9   0.0   0.0   0.0   0.0   0.0   0.3   0.8   0.9   1.0   1.0
    ##  [8,]   0.2   0.0   0.0   0.0   0.0   0.0   0.0   0.3   0.7   0.8   0.8
    ##  [9,]   0.1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.2
    ##       [,91] [,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100]
    ##  [1,]   0.7   1.0   1.0   1.0   1.0   1.0   1.0   1.0   1.0      1
    ##  [2,]   0.3   0.9   1.0   1.0   1.0   1.0   1.0   1.0   1.0      1
    ##  [3,]   0.0   0.0   0.7   1.0   1.0   1.0   1.0   1.0   1.0      1
    ##  [4,]   0.0   0.1   0.2   0.4   0.9   0.9   0.9   0.9   1.0      1
    ##  [5,]   0.0   0.0   0.1   0.6   1.0   1.0   1.0   1.0   1.0      1
    ##  [6,]   0.0   0.0   0.0   0.0   0.0   0.7   0.9   1.0   1.0      1
    ##  [7,]   0.0   0.0   0.0   0.0   0.0   0.3   0.8   0.9   1.0      1
    ##  [8,]   0.0   0.0   0.0   0.0   0.1   0.1   0.4   0.9   1.0      1
    ##  [9,]   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.2   0.8      1

Estimate the true signal
------------------------

### Construct bootstrap samples

Sample from the columns of `R.input` with replacement in order to construct `num.boot` bootstrap rank matrices.

``` r
num.boot = 50  # number of bootstrap matrices
boots = list()  # list of bootstrap matrices
boots[[1]] = R.input
for (b in 2:(num.boot + 1)) {
    set.seed(b)
    ind = sample.int(n, n, replace = TRUE)
    boots[[b]] = R.input[, ind]
}
```

### Run the Metropolis MCMC algorithm on each bootstrap matrix.

**!WARNING!** This step takes longer time, depending on the number of cores used. Maximum efficiacy on a single machine can be achieved by using `num.sim` cores, as the individual chains are run in parallel. The results of this step are saved in `Example_all_results.Rdata`.

``` r
chain.length = 10000 # length of an MCMC chain
num.sim = 10 # number of independent chains
cores = detectCores() # number of cores to run the algorithm on
res = list() # here the results are stored

set.seed(123)
theta.inis = vector('list',num.sim)# initial guesses - one for each chain
theta.inis = lapply(theta.inis, function(x) x= runif(p, -1,1))
runtime = numeric()
for (b in 1:(num.boot+1)) ## for each bootstrap matrix
{
  runtime[b] = system.time({
  res.temp = MCMC.metropolis(theta.inis = theta.inis, in.data=list(boots[[b]], F.input), dev=0.1, chain.len=chain.length, cores=cores)
  })[3]
  res[[b]] = res.temp
  save(runtime, res.temp, file = paste0('Example_boot',b,'.RData'))
}
save(res, file = 'Example_all_results.Rdata')
```

Gather the results
------------------

If the previous step was performed, comment out the first two lines.

``` r
load("Example_all_results.Rdata")
num.sim = 10
indiv.estimates = replicate(num.boot + 1, matrix(nrow = p, ncol = num.sim), 
    simplify = FALSE)  # initiate list of matrices, where the solutions are stored
main.nrm.estimates = matrix(nrow = p, ncol = num.boot + 1)

for (b in 1:(num.boot + 1)) {
    for (i in 1:num.sim) {
        if (is.numeric(res[[b]]$avg[[i]]$x.in.min)) 
            indiv.estimates[[b]][, i] = res[[b]]$avg[[i]]$x.in.min else indiv.estimates[[b]][, i] = res[[b]]$avg[[i]]$x.in.min[, 
            1]
    }
    m = apply(indiv.estimates[[b]], 1, median)  # median over 10 independent chains
    m.nrm = m/norm_vec(m)  # normalisation
    
    main.nrm.estimates[, b] = m.nrm
}
```

Plot the results with ggplot
----------------------------

### Prepare

``` r
melt.data = data.frame(melt(main.nrm.estimates))
colnames(melt.data) = c("object", "boot", "value")

truth.est = data.frame(estimate = main.nrm.estimates[, 1], true.signal = theta.true.norm, 
    object = 1:nrow(main.nrm.estimates))
truth.est = data.frame(melt(truth.est, id.vars = "object"))
```

### Plot

``` r
pl = ggplot(melt.data, aes(factor(object), value)) + ylim(-1, 
    1)
pl = pl + xlab("Object") + ylab("Signal value") + geom_violin(trim = TRUE, 
    scale = "count", size = 0.5)  ## violin plots
pl = pl + geom_line(data = truth.est, aes(object, value, col = variable), 
    size = 1)  # the true and estimated signal
pl = pl + stat_summary(data = melt.data, fun.data = mean_sterr, 
    geom = "errorbar", color = "green", size = 0.2)  # error bars
pl = pl + stat_summary(fun.y = mean, geom = "point", fill = "black", 
    shape = 21, size = 2, position = position_dodge(width = 0.9))
pl = pl + theme_bw() + theme(text = element_text(size = 15)) + 
    theme(legend.title = element_blank())  # white background, big letters
pl = pl + scale_colour_manual(values = c("blue", "red"), labels = c("Estimate", 
    "True signal"))  # colors for true and estimated
pl = pl + theme(legend.text = element_text(size = 15), legend.key = element_blank(), 
    legend.position = c(0.75, 0.65), legend.justification = c(0, 
        0), legend.margin = unit(0, "cm"), panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank())
pl + scale_x_discrete(breaks = 1:p, labels = paste0("o", 1:p))  # customised ticks
```

![](Example_markdown_files/figure-markdown_github/unnamed-chunk-9-1.png)

#############################################################################
#############################################################################
#######  The MultiRankS ALGORITHM
#######  from the paper by Svendova, Schimek 2017:
#######  A novel method for estimating the common signals for consensus across multiple ranked lists
#############################################################################

# Functions' description:
#    norm_vec: normalises a vector in order to have unit norm
#    st.error.estim: calculates standard error
#    mean_sterr: calculates standard error and the mean of the values
#    F_perm: calculates rank probability matrix for window height of specified maximum
#    J_theta: evaluates the objective function for a specified probability matrix
#    FRgeneral: given a vector and input matrix, calculates a rank matrix most similar to the input one
#    what.is.suitable.sigma: estimates a standard deviation needed to create the most similar rank matrix (see FRgeneral)
#    MCMC.metropolis: runs onerun() function as many times as there are specified initial guesses (theta.inis). Distributes the run into a specified number of cores.
#    onerun: runs the function run_metropolis_MCMC() and summarises its results
#   run_metropolis_MCMC: runs classical Metropolis-Hastings algorithm with a specified fixed step size
#    proposalfunction: calculates a next step proposal in Metropolis-Hastings algorithm
#    adaptive.MCMC.metropolis: runs onerun.adaptive() function as many times as there are specified initial guesses (theta.inis). Distributes the run into a specified number of cores.
#    onerun.adaptive: runs the function MCMC() of package adaptMCMC, i.e. adaptive Metropolis algorithm, as described in Vihola, M. (2012). Robust adaptive metropolis algorithm with coerced acceptance rate. Statistics and Computing, 22(5):997–1008.

library(gtools)
library(reshape)
library(ggplot2)
library(grid)
library(parallel)
library(adaptMCMC)



gather.results=function(res)
{
  num.boot = length(res) - 1
  num.sim = length(res[[1]]$avg)
  p = length(res[[1]]$avg[[1]]$x.in.min)
  indiv.estimates = replicate(num.boot+1, matrix(nrow=p, ncol=num.sim), simplify = FALSE) # initiate list of matrices, where the solutions are stored
  main.nrm.estimates = matrix(nrow=p, ncol=num.boot+1)
  
  for (b in 1:(num.boot+1))
  {
    for (i in 1:num.sim)
    {
      if (is.numeric(res[[b]]$avg[[i]]$x.in.min))
        indiv.estimates[[b]][,i] = res[[b]]$avg[[i]]$x.in.min
      else indiv.estimates[[b]][,i] = res[[b]]$avg[[i]]$x.in.min[,1]
    }
    m = apply(indiv.estimates[[b]], 1, median) # median over 10 independent chains
    m.nrm = m/norm_vec(m) # normalisation
    
    main.nrm.estimates[,b] = m.nrm
  }
  est.plus.SE = data.frame(signal.estimate = apply(main.nrm.estimates,1,mean), SE=apply(main.nrm.estimates,1,st.error.estim))
  return(est.plus.SE)
}

run.adaptiveMCMC.example = function(F.input, boots, num.sim, chain.length)
{
  num.boot = length(boots)-1
  p=nrow(boots[[1]]); n=ncol(boots[[1]])
  cores = detectCores() # number of cores to run the algorithm on
  res = list() # here the results are stored
  
  set.seed(123)
  theta.inis = vector('list',num.sim)# initial guesses - one for each chain
  theta.inis = lapply(theta.inis, function(x) x= runif(p, -1,1))
  for (b in 1:(num.boot+1)) ## for each bootstrap matrix
  {
      res.temp = adaptive.MCMC.metropolis(theta.inis = theta.inis, in.data=list(boots[[b]], F.input), chain.len=chain.length, cores=cores, l_max=length(F.input))
    res[[b]] = res.temp
  }
  return(res)
}



generate.bootstrap.samples = function(R,num.boot)
{
  # R - the input matrix
  # num.boot - number of bootstrap matrices
  n = ncol(R); p=nrow(R)
  boots = list() # list of bootstrap matrices
  boots[[1]] = R 
  for (b in 2:(num.boot+1))
  {
    set.seed(b)
    ind = sample.int(n, n, replace=TRUE)
    boots[[b]] = R[,ind]
  }
  return(boots)
}

generate.random.rank.matrix = function(p,n)
{
  set.seed(123)
  sigmas = abs(rnorm(n,0,0.4)) # random standard deviations for each assessors 
  X.input = matrix(nrow=p, ncol=n) # the matrix of observed attributes
  theta.true = rnorm(p,0,1) # true signal
  for (i in 1:n)
  {
    set.seed(i)
    X.input[,i] = theta.true + rnorm(p,0,sigmas[i])
  }
  
  R.input = matrix(nrow=p, ncol=n) # the observed rankings 
  R.input = apply(X.input, 2, function(x) rank(-x)) 
  return(list(R=R.input,theta=theta.true))
}

norm_vec = function(x) sqrt(sum(x^2)) # vector normalisation

st.error.estim = function(x){sqrt(sum((x-mean(x))^2/(length(x)-1)))} # standard error estimation

mean_sterr <- function(x) { # estimation of mean and 2SE from numeric vector x
  data.frame("y" = mean(x), "ymin" = mean(x) - 2*st.error.estim(x), "ymax" = mean(x) + 2*st.error.estim(x))
}


F_perm <- function(R, l_max){
  ## Calculates a probability matrix for each window size up to l_max. Outputs a list of all the probability matrices.
  
  ## R - matrix of rankings, objects in rows, assessors in columns
  ## l_max - maximum window height
  
  p=nrow(R) # number of objects
  n=ncol(R) # number of assessors
  
  F.l = list() ## list of matrices F, each for different l
  for (ell in 1:l_max)
  { 
    s_mat <- unname(as.matrix(do.call(expand.grid,rep(list(seq_len(p)),ell))))
    F_perm <- matrix(NA_real_,p-ell+1,nrow(s_mat))
    for (sri in seq_len(nrow(s_mat))) {
      s <- s_mat[sri,];
      F_perm[,sri] <- rowSums(Reduce(`&`,Map(function(e,i) R[i:(p-length(s)+i),]<=e,s,seq_along(s))))/n 
    }
    
    F.l[[ell]] = F_perm
    names(F.l)[[ell]] = paste("l=",ell,sep="") 
  }
  return(F.l = F.l)
}

J_theta = function(Fi, F.temp, l_max)
{
  ## Calulates the value of the objective function J.
  
  # Fi - the input list of probability matrices (calculated from the true rank matrix)
  # F.temp - another list of probability matrices (calculated from a current estimate)
  SumSl = 0
  Suml = 0
  for (ell in 1:l_max)
  {
    for (sl in 1:ncol(Fi[[ell]]))
    {
      F1 = Fi[[ell]][,sl]
      F2 = F.temp[[ell]][,sl] 
      SumSl = SumSl + (F1 - F2)%*%(F1 - F2)
    }
    Suml = Suml + SumSl
  }
  Suml = Suml / l_max ## normalization 
  
  return(Suml)
}


FRgeneral <- function(theta, studies, l_max, R.input, F.input,increments=NULL) 
{
  ## Calculates a rank matrix R.temp, based on numerical vector theta. Matrix R.temp should resemble the input rank matrix R.input. The resemblance is achieved by generating a 'suitable' st.dev. with function what.is.suitable.sigma(). Outputs the rank matrix R.temp, its probablity matrices F.temp and the value of J J.theta.
  
  # theta - a numerical vector (current estimate)
  # l_max - maximum window size
  # R.input - input rank matrix
  # F.input - input list of probability matrices of R.input
  # increments - see what.is.suitable.sigma() function
  
  sigma = what.is.suitable.sigma(theta,R.input,increments)
  X = matrix(nrow=length(theta), ncol=studies)
  X = apply(X, 2, function(x) theta + rnorm(nrow(R.input), mean=0, sigma))
  R.temp = apply(X, 2, function(x) rank(-x, ties.method='random'))
  
  F.temp = F_perm(R.temp,l_max)
  J.theta = J_theta(F.input, F.temp, l_max)
  return(list(J = J.theta, F = F.temp, R = R.temp))	
}

what.is.suitable.sigma <- function(theta, R.input, increments=NULL)
{
  ## Estimates a st.dev. needed to create a rank matrix as similar to R.input as possible, using the vector theta. Outputs the most suitable sigma for theta.
  
  # theta - a numerical vector (current estimate)
  # R.input - input rank matrix
  # increments - the size of the increments of sigmas tested for suitability. Larger increments makes the algorithm faster, but tests less sigma-posibilities. Should be adjusted depending on the number of objects p (for p~10 increments = 0.05 are satisfactory)
  
  p = nrow(R.input); n = ncol(R.input)
  if (is.null(increments)) increments = 0.05
  corInput = cor(R.input, method='spearman')
  rho = median(corInput[lower.tri(corInput)]) ## average correlation
  sigma_range = seq(0.01,1, by = increments) ## testing what correlation will be caused by each of the sigmas
  rho.temp=numeric()
  Xs = list() # list of measurement matrices, as produces by using the individual sigmas from sigma_range
  for (sig in 1:length(sigma_range))
  {
    Xs[[sig]] = apply(matrix(nrow=p, ncol=n),2,function(x) x=theta + rnorm(p, mean=0, sd=sigma_range[sig])) 
  }
  corm=Map(function(x) cor(x,method='spearman'), Xs)  
  rho.temp = sapply(corm, function(x) median(x[lower.tri(x)]))
  sigma = sigma_range[which(abs(rho.temp-rho)==min(abs(rho.temp-rho)))][1] 
  return(sigma)
}

real.fitness <- function(theta, n, p, l_max, F.input, R.input, increments=NULL) ### DEFINE FITNESS/OBJECTIVE FUNCTION
{
  if (is.null(increments)) increments=0.01
  ### calculating appropriate sigma for each iteration:
  sigma = what.is.suitable.sigma(theta,R.input,increments)
  Z_ij = matrix(nrow=p, ncol=n)
  Z_ij = apply(Z_ij, 2, function(x) rnorm(p, mean=0, sigma))
  X = apply(Z_ij, 2, function(x) theta + x)
  R = apply(X, 2, function(x) rank(-x, ties.method='random'))
  F.temp = F_perm(R,l_max)
  J.theta = J_theta(F.input, F.temp, l_max)
  return(-J.theta)
} 

#-------------------------------------------
#
#        Classical Metropolis algorithm 
#
#-------------------------------------------

MCMC.metropolis <- function(theta.inis,in.data, dev, chain.len, cores=detectCores(), increments=NULL){
  ## Runs several independent Metropolis MCMC chains (function onerun()). Number of chains depends on the length of initial guesses theta.inis. Outputs list of the results (chain + best estimate + minimum found) and runtime
  
  # theta.inis - a list of vectors of initial guesses, the length of the list determines the number of chains 
  # in.data - list of 2: first rank matrix R, second list of matrices F (result of F_perm)
  # dev - the standard deviation defining every next proposal step of the random walk
  # chain.len - length of each chain
  # cores - number of cores to be used
  # increments - see what.is.suitable.sigma() function
  
  # clusterApply() for Windows
  if (Sys.info()[1] == "Windows"){
    cl <- makeCluster(cores)
    clusterExport(cl, list("onerun", "run_metropolis_MCMC", "FRgeneral", "what.is.suitable.sigma", "F_perm", "J_theta", "proposalfunction"))
    runtime <- system.time({
      res = clusterApplyLB(cl=cl, x=theta.inis, fun=onerun, input=in.data, dev=dev, its=chain.len, l_max=2)
    })[3]
    stopCluster(cl) # Don't forget to do this!
    
    # mclapply() for everybody else
  } else {
    runtime <- system.time({
      res = mclapply(X=theta.inis, FUN=function(x) onerun(theta.ini=x, input=in.data, dev=dev, its=chain.len, l_max=2, increments), mc.cores=cores)
    })[3]
  }
  return(list(avg=res, runtime=runtime))
}

onerun <- function(theta.ini, input, dev=NULL, its, l_max, increments=NULL){ # 
  ## Runs one chain of the Metropolis MCMC (function run_metropolis_MCMC()) and finds the minima. Outputs all found vectors where the minimum was found (x.in.min), all points of the chain(res.chain), and the value of the minimum (minJ).
  
  # theta.ini - the initial guess (numeric vector)
  # input - list of 2: first rank matrix R, second list of matrices F (result of F_perm)
  # dev - the standard deviation defining every next proposal step of the random walk
  # its - length of the chain 
  # increments - see what.is.suitable.sigma() function
  if(is.null(dev)) dev=0.1
  if (is.null(increments)) increments=0.01
  if (is.null(l_max)) l_max=2
  
  res.chain = list()
  mins = list()
  x.at.mins = list()
  acceptance = list()
  
  p = nrow(input[[1]])
  Ri = input[[1]]
  Fi = input[[2]]
  
  res.chain = run_metropolis_MCMC(theta.ini, iterations=its, dev=dev, Ri=Ri, Fi=Fi, l_0=l_max, increments)
  acceptance = 1-mean(duplicated(t(res.chain$chain))) 
  if (acceptance < 0.15) warning(paste0('The acceptance rate is below 15% (',acceptance,'). The step size probably needs to be smaller.'))
  if (acceptance > 0.5) warning(paste0('The acceptance rate is above 50% (',acceptance,'). The step size probably needs to be larger.'))
  mins = which(res.chain$Js == min(res.chain$Js))
  x.at.mins = res.chain$chain[,mins]
  if (is.matrix(x.at.mins)) 
  { all.x.in.mins = as.data.frame(x.at.mins)
  dupl = apply(all.x.in.mins,1,duplicated) ## find duplicated columns
  all.x.in.mins = all.x.in.mins[,which(apply(dupl,1,function(x) sum(x)) !=p)] # keep only unique columns
  }
  else all.x.in.mins = x.at.mins
  return(list(x.in.min=all.x.in.mins, chain=res.chain, minJ=min(res.chain$Js), acceptance=acceptance))
}

run_metropolis_MCMC <- function(theta.ini, iterations, dev, Ri, Fi,l_0,increments=NULL){
  ## Calculates Metropolis MCMC chain. Outputs the chain (all visited points) and the values of J in ech point.
  
  # theta.ini - initial point (numerical vector)
  # iterations - number of mcmc steps (length of the chain)
  # dev - the standard deviation defining every next proposal step of the random walk
  # Ri - input rank matrix
  # Fi - input list of probability matrices of Ri
  # increments - see what.is.suitable.sigma() function
  if (is.null(dev)) dev=0.1
  if (is.null(increments)) increments=0.01
  if (is.null(l_0)) l_0=2
  
  p=nrow(Ri); n=ncol(Ri)
  chain = array(dim = c(length(theta.ini),iterations+1)) # array of points I visit
  F.for.chain = list() # saving the F values (might be needed for some reason)
  Js = numeric() # J for each point I visit
  chain[,1] = theta.ini
  FR_ini = FRgeneral(theta.ini, n, l_0, Ri, Fi, increments) # calculates the F and J for the initial point
  F.for.chain[[1]] = FR_ini$F
  Js[1] = FR_ini$J
  
  for (i in 1:iterations){
    proposal = proposalfunction(param=chain[,i], stdev=dev) # create next point
    FR_proposal = FRgeneral(theta=proposal, studies=n, l_max=l_0, R.input=Ri, F.input=Fi,increments) # calculates the F and J for a current point
    J_proposal = FR_proposal$J
    
    probab = min(1,exp(-(J_proposal - Js[i]))) # the probability of acceptance 
    
    if (runif(1) < probab){ ## accept the new point with probability probab
      chain[,i+1] = proposal
      Js[i+1] = J_proposal
      F.for.chain[[i+1]] =  FR_proposal$F
    }else{ ## otherwise stay in the current point
      chain[,i+1] = chain[,i]
      Js[i+1] = Js[i]
      F.for.chain[[i+1]] =  F.for.chain[[i]]
    }
  }
  return(list(chain=chain, Js=Js))
}

proposalfunction <- function(param, stdev){
  ## Calculates a proposal of the next step in MCMC algorithm
  # param - a numerical vector (current point)
  # stdev - st.dev. of the random values added the the previous point to create the next point 
  
  e = rnorm(length(param),mean = 0, sd = stdev)
  return(param+e)
}



#-------------------------------------------
#
#        Adaptive Metropolis algorithm 
#        with usage of the package adaptMCMC
#
#-------------------------------------------

adaptive.MCMC.metropolis <- function(theta.inis,in.data, chain.len,  cores=detectCores(), l_max=NULL, increments=NULL, scale=NULL){
  ## Runs several independent Metropolis MCMC chains (function onerun()). Number of chains depends on the length of initial guesses theta.inis. Outputs list of the results (chain + best estimate + minimum found) and runtime
  
  # theta.inis - a list of vectors of initial guesses, the length of the list determines the number of chains 
  # in.data - list of 2: first rank matrix R, second list of matrices F (result of F_perm)
  # dev - the standard deviation defining every next proposal step of the random walk
  # chain.len - length of each chain
  # cores - number of cores to be used
  # increments - see what.is.suitable.sigma() function
  if (is.null(increments)) increments=0.01
  if (is.null(l_max)) l_max=2
  if (is.null(scale)) scale=rep(0.1,length(theta.inis[[1]]))
  # clusterApply() for Windows
  if (Sys.info()[1] == "Windows"){
    cl <- makeCluster(cores)
    clusterExport(cl, list("onerun.adaptive", "FRgeneral", "what.is.suitable.sigma", "F_perm", "J_theta", "real.fitness"))
    runtime <- system.time({
      res = clusterApplyLB(cl=cl, x=theta.inis, fun=onerun.adaptive, input=in.data, its=chain.len, l_max=l_max, scale=scale, increments=increments)
    })[3]
    stopCluster(cl) # Don't forget to do this!
    
    # mclapply() for everybody else
  } else {
    runtime <- system.time({
      res = mclapply(X=theta.inis, FUN=function(x) onerun.adaptive(theta.ini=x, input=in.data, its=chain.len, l_max=l_max, scale=scale, increments=increments), mc.cores=cores)
    })[3]
  }
  return(list(avg=res, runtime=runtime))
}

onerun.adaptive <- function(theta.ini, input,  its, l_max=NULL, increments=NULL, scale=NULL){ # 
  ## Runs one chain of the Metropolis MCMC (function run_metropolis_MCMC()) and finds the minima. Outputs all found vectors where the minimum was found (x.in.min), all points of the chain(res.chain), and the value of the minimum (minJ).
  
  # theta.ini - the initial guess (numeric vector)
  # input - list of 2: first rank matrix R, second list of matrices F (result of F_perm)
  # its - length of the chain
  # l_max - maximum height of the sliding window
  # increments - see what.is.suitable.sigma() function
  # scale - starting step size before adapting
  if (is.null(increments)) increments=0.01
  if (is.null(l_max)) l_max=2
  if (is.null(scale)) scale=rep(0.1,length(theta.ini))
  
  res.chain = list()
  mins = list()
  x.at.mins = list()
  acceptance = list()
  
  p = nrow(input[[1]])
  Ri = input[[1]]
  Fi = input[[2]]
  
  res = MCMC(p=function(x) real.fitness(x,n=ncol(Ri),p=nrow(Ri),l_max,Fi,Ri,increments), n=its, init=theta.ini, acc.rate = 0.234, scale=scale,adapt=TRUE, n.start=100)
  minimum = res$samples[its,]
  minJ=-max(res$log.p)
  return(list(x.in.min=minimum, Js=res$log.p, minJ=minJ))
}

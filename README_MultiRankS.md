## MultiRankS_funs.R 

 Coded by Vendula Svendova


Functions description:

  - norm_vec: normalises a vector in order to have unit norm
  - st.error.estim: calculates standard error
  - mean_sterr: calculates standard error and the mean of the values
  - F_perm: calculates rank probability matrix for window height of specified maximum
  - J_theta: evaluates the objective function for a specified probability matrix
  - FRgeneral: given a vector and input matrix, calculates a rank matrix most similar to the input one
  - what.is.suitable.sigma: estimates a standard deviation needed to create the most similar rank matrix (see FRgeneral)
  - MCMC.metropolis: runs onerun() function as many times as there are specified initial guesses (theta.inis). Distributes the run into a specified number of cores.
  - onerun: runs the function run_metropolis_MCMC() and summarises its results
  - run_metropolis_MCMC: runs classical Metropolis-Hastings algorithm with a specified fixed step size
  - proposalfunction: calculates a next step proposal in Metropolis-Hastings algorithm
  - adaptive.MCMC.metropolis: runs onerun.adaptive() function as many times as there are specified initial guesses (theta.inis). Distributes the run into a specified number of cores.
  - onerun.adaptive: runs the function MCMC() of package adaptMCMC, i.e. adaptive Metropolis algorithm, as described in Vihola, M. (2012). Robust adaptive metropolis algorithm with coerced acceptance rate. Statistics and Computing, 22(5):997–1008.
  

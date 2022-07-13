get_values <- function(D, y, z, lambda){
  ### Computes objective function by each term (stress, regularization) given a configuration
  dist_z <- dist(z)
  dist_z <- as.matrix(dist_z)
  stress <- 0.5*sum((D-dist_z)^2)
  stress.between <- sum(((D - dist_z)[y==1,y==2])^2)
  stress.within <- 0.5*(sum(((D - dist_z)[y==1,y==1])^2) + sum(((D - dist_z)[y==2,y==2])^2))
  return(list(stress=, stress.within=, stress.between=, regularization=, ))
}

noname_once <- function(D=NULL, y, lambda=.5, m=2, x=NULL){
  ### Performs gradient descent to minimize objective function (stress plus regularization)
  return(list(A=, B=, ))
}

noname <- function(D=NULL, y, lambda=.5, m=2, x=NULL, n_rep=5){
  ### Iterates noname_once with randomly starting z
  # Inputs
  # D: distance matrix, n by n
  # y: binary class labels, Either 0 or 1
  # lambda: regularization parameter. Positive real number
  # m: target number of reduced dimension. Integer less than n
  # x: design matrix, n by p
  # Outputs
  # z: optimized configuration. Matrix n by m
  return()
}
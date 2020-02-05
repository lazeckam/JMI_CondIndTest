library(mvnfast)
library(infotheo)
library(e1071)

# Datasets generation ----

generate_GTM <- function(k=10, n=1000, gamma=2/3, Sigma=diag(k), nbins=2){
  
  No_ones <- rbinom(1, n, prob=1/2)
  
  X <- rbind(rmvn(No_ones, gamma^(0:(k-1)), Sigma),
             rmvn(n - No_ones, rep(0, k), Sigma))
  
  X <- data.frame(apply(X, 2, function(x) infotheo::discretize(x, nbins=nbins))) - 1
  
  X <- data.frame(X, as.vector(as.matrix(infotheo::discretize(rnorm(n, X[ , 1], 1), nbins=nbins))) - 1)
  colnames(X) <- c(paste0("X", 1:k)," X1(1)")
  
  list(Y=c(rep(1, No_ones), rep(0, n - No_ones)), X=data.frame(X))
}


generate_logist_0 <- function(k=10, n=1000, nbins=3){
  
  X <- matrix(rnorm(n*k, 0, 1), n, k)
  
  # 2 bins
  X <- data.frame(apply(X, 2, function(x) 2*infotheo::discretize(x, nbins=2))) - 3
  X <- data.frame(X, as.vector(as.matrix(2*infotheo::discretize(rnorm(n, X[ , 1], 1), nbins=2))) - 3)
  
  # 3 bins
  # X <- data.frame(apply(X, 2, function(x) infotheo::discretize(x, nbins=3))) - 2
  # X <- data.frame(X, as.vector(as.matrix(infotheo::discretize(rnorm(n, X[ , 1], 1), nbins=3))) - 2)
  
  colnames(X) <- c(paste0("X", 1:k)," X1(1)")
  
  Y <- 1*(runif(n, 0, 1) <  exp(apply(data.frame(X[,1:k]), 1, sum))/(1 + exp(apply(data.frame(X[,1:k]), 1, sum))))
  
  # 2 bins
  X <- data.frame(apply(X, 2, function(x) ifelse(x==-1,0,1)))
  
  # 3 bins
  # X <- data.frame(apply(X, 2, function(x) ifelse(x==-1,2,1)))
  
  colnames(X) <- c(paste0("X", 1:k)," X1(1)")
  
  list(Y=Y, X=data.frame(X))
}

generate_logist_1 <- function(k=10, n=1000, gamma=1){
  
  X <- matrix(rnorm(n*(k+1), 0, 1), n, k+1)
  
  # 2 bins
  X <- data.frame(apply(X, 2, function(x) 2*infotheo::discretize(x, nbins=2))) - 3
  
  # 3 bins
  # X <- data.frame(apply(X, 2, function(x) infotheo::discretize(x, nbins=3))) - 2
  
  p <- exp(apply(data.frame(X[,1:k]), 1, sum) + gamma*X[,1]*X[,k+1])
  Y <- 1*(runif(n, 0, 1) <  p/(1 + p))
  
  # 2 bins
  X <- data.frame(apply(X, 2, function(x) ifelse(x==-1,0,x)))
  
  # 3 bins
  # X <- data.frame(apply(X, 2, function(x) ifelse(x==-1,2,x)))
  
  colnames(X) <- c(paste0("X", 1:k)," X1(1)")
  
  list(Y=Y, X=data.frame(X))
}


# Permutation over layers ----

permute <- function(X, Y, Z){
  
  Z <- as.matrix(Z)
  
  if (ncol(Z) > 1) {
    Z_one_column <- apply(Z , 1, paste, collapse = ",")
  } else {
    Z_one_column <- Z
  }
  
  Z_layers <- unique(Z_one_column)
  X_perm <- numeric(length(X))
  
  for(layer in Z_layers){
    temp <- X[Z_one_column == layer]
    temp <- temp[sample(1:length(temp), length(temp), replace=FALSE)]
    X_perm[Z_one_column == layer] <- temp
  }
  
  return(list(X=X_perm, Y=Y, Z=Z))
}

# Other ----

covariance_matrix <- function(n, rho){
  
  H <- abs(outer(1:n, 1:n, "-")) 
  sigma <- rho^H
  
  return(sigma)
}

zhang_params <- function(x){
  
  mu <- mean(x)
  sigma2 <- var(x)
  mu3 <- moment(x,order = 3, center = TRUE)
  
  if(sigma2 > 0 & mu3!=0){
    alpha <- 1/4 * mu3/sigma2
    d <- 8 * sigma2^3/mu3^2
    beta <- mu - alpha*d
  }else{
    if(mu>0){
      alpha <- 1
      d <- mu
      beta <- 0
    }else if(mu<0){
      alpha <- -1
      d <- -mu
      beta <- 0
    }else{
      alpha <- 1
      d <- 1
      beta <- 0
    }
  }
  return(list(alpha = alpha, d = d, beta = beta))
}
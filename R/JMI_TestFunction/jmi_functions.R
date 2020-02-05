require(e1071)

# JMI Conditional Independence Test: ----

permute_sample <- function(X, Z_one_column){
  
  # Output: a vector X permuted on Z_one_column layers
  
  Z_layers <- unique(Z_one_column)
  X_perm <- numeric(length(X))
  
  for(layer in Z_layers){
    temp <- X[Z_one_column == layer]
    temp <- temp[sample(1:length(temp), length(temp), replace=FALSE)]
    X_perm[Z_one_column == layer] <- temp
  }
  
  return(X_perm)
}

JMI_components <- function(X, Y, Z){
  
  # Output: a vector containing values of sum_j log p(x_i, y_i, z_ij)*p(z_ij)/p(x_i, z_ij)*p(y_i, z_ij)
  #         i - an index of observations, j - an index of conditioning variables

  X <- as.numeric(as.factor(as.numeric(X)))
  Y <- as.numeric(as.factor(as.numeric(Y)))
  Z <- data.frame(Z)
  n <- length(X)

  z <- apply(Z, 2, as.character)
  
  w <- rep(0, n)
  
  for(j in 1:ncol(Z)){
    tab_xy_cond_z <- table(X, Y, factor(Z[,j]))
    tab_xy_cond_z <- apply(tab_xy_cond_z, 3, function(t) {
      p <- t/sum(t)
      p_indep <- apply(p, 1, sum) %o% apply(p, 2, sum)
      data.frame(log0(as.matrix(p/p_indep)))
      })
    for(i in 1:n){
      w[i] <- w[i] + tab_xy_cond_z[[z[i,j]]][X[i], Y[i]]
    }
    
  }
  
  return(w)
}

zhang_params <- function(x){
  
  mu <- mean(x)
  sigma2 <- var(x)
  mu3 <- moment(x, order = 3, center = TRUE)
  
  if(sigma2 > 0 & mu3!=0){
    alpha <- 1/4 * mu3/sigma2
    d <- 8 * sigma2^3/mu3^2
    beta <- mu - alpha*d
  }else{
    if(mu > 0){
      alpha <- 1
      d <- mu
      beta <- 0
    }else if(mu < 0){
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

log0 <- function(p){
  
  p <- as.matrix(p)
  
  res <- log(p)
  res[!is.finite(res)] <- 0
  
  return(res)  
}

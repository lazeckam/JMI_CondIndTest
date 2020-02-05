
##### Packages & functions ----

# Computation of:
# JMI_component_i = sum_i log p(x, y, z_i)*p(z_i)/p(x, z_i)*p(y, z_i), i - observation index
# JMI = mean(JMI_component_i)
library(Rcpp)
sourceCpp("JMI_components.cpp")

# Generating datasets functions & others:
source("functions.R")

##### Grid of parameters ----

params <- list(gamma = c(1),         # parameter that controls how much information X_i contains about Y
               rho = c(0),           # strength of dependence in X
               k = 2:20,             # number of conditioning variables
               n = c(500, 1000),     # number of observations
               dep = c(TRUE, FALSE)) # dependence
grid.params1 <- expand.grid(params)

params <- list(gamma = c(0.5, 1),
               rho = seq(0.8, 0, -0.2),
               k = 4,
               n = c(1000),
               dep = c(TRUE))
grid.params2 <- expand.grid(params)

grid.params <- rbind(grid.params1, grid.params2)
n.params    <- nrow(grid.params)
all.params  <- colnames(grid.params)

##### Other parameters ----

N <- 500   # number of experiments
B <- 50    # number of permutations

choose.distribution <- rep(0, N)     # p-values of the "swich" test
p.val.distribution  <- rep(0, N)     # p-values of norm/chi^2_scale tests
                                     # alpha is set as 0.05
results <- vector("list", n.params) 

##### Main loop ----

for(i in 1:n.params){
  
  p.val.distribution.norm <- rep(NA, N)
  p.val.distribution.chi <- rep(NA, N)
  p.val.distribution.chi.scale <- rep(NA, N)
  
  for(par.name.tmp in all.params) assign(par.name.tmp, grid.params[i, par.name.tmp])
  
  Sigma <- covariance_matrix(k, rho)
  
  cat(paste0(c("|",rep("_", 50), "|", i, "/", n.params), collapse=""), "\n")
  cat("|")
  gwiazdki <- 0
  
  for(j in 1:N){
    
    # Generating data from the GTM (Generative Tree Model)
    #      X indep. Y | Z  will be tested
    
    df <- generate_GTM(k=k, n=n, gamma=gamma, Sigma=Sigma)
    Y <- df$Y
    Z <- as.matrix(df$X[, 2:k])
    
    X <- df$X[, k+1]     # (H_{02}) dep == TRUE  -> X = X_1^{(1)} (M1)
    if(!dep){            # (H_{03}) dep == FALSE -> X = X_1       (M1)
      X <- df$X[, 1] 
    }
    
    # JMI variance computations
    
    var.values <- JMI_components(X, Y, Z)  # JMI components values
    mean.value <- mean(var.values)         # JMI value (mean of components)
    
    var.mean.full <- replicate(B, {        # variance & mean values for the full dataset with permuted X over Z=z layers
      X_perm <- permute(X, Y, Z)$X
      var_values_tmp <- JMI_components(X_perm, Y, Z)
      c(var(var_values_tmp), mean(var_values_tmp))
    })
    
    variance.distribution <- var.mean.full[1,]
    mean.distribution <- var.mean.full[2,] # JMI values for full set with X indep. Y|Z (H_0)
    var.full <- var.mean.full[1,]          # JMI variances for full set with X indep. Y|Z (H_0)
    
    var.half <- replicate(B, {             # variance values for half of the dataset with permuted X over Z=z layers
      X_perm <- permute(X, Y, Z)$X
      i <- sample(1:n, n/2)                # random half of the dataset, different for every loop
      var(JMI_components(X_perm[i], Y[i], as.matrix(Z[i,])))
    })
    
    test <- t.test(var.full, var.half, "less")  # "swich" test: variance comparison
    choose.distribution[j] <- test$p.value      # H_0 ~ normal, !H_0 ~ chi^2
    
    if(choose.distribution[j] > 0.05 | mean(2*n*mean.distribution) <= 0){         
                                           # H_0 of the "swich" test
                                           #  or exception to make sure that the mean in chi^2 is greater than 0
      
      p.val.distribution.norm[j] <- pnorm(sqrt(n)*mean.value, 
                                          mean=sqrt(n)*mean(mean.distribution), 
                                          sd=sqrt(mean(variance.distribution)), lower.tail=FALSE)
      
      p.val.distribution[j] <- p.val.distribution.norm[j]
      
    }else{
      
      p.val.distribution.chi[j] <- pchisq(2*n*mean.value, 
                                          df=mean(2*n*mean.distribution), lower.tail = FALSE)
      
      par <-  zhang_params(2*n*mean.distribution)
      lambda <- par$d
      alpha <- par$alpha
      beta <- par$beta
      p.val.distribution.chi.scale[j] <- pchisq((2*n*mean.value - beta) / alpha, 
                                                df=lambda, lower.tail = FALSE)
      
      p.val.distribution[j] <- p.val.distribution.chi.scale[j]
    }
    
    if(floor(50*j/N - gwiazdki) >= 1){
      replicate(floor(50*j/N - gwiazdki), cat("*"))
      gwiazdki <- gwiazdki+floor(50*j/N - gwiazdki)
    }
    
  }
  cat("|\n")
  
  results[[i]] <- list(parameters = grid.params[i,],                          
                       p.val.choose.distribution = choose.distribution,
                       choose.distribution.norm = choose.distribution > 0.05,
                       p.val.distribution = p.val.distribution,
                       p.val.distribution.norm = p.val.distribution.norm,
                       p.val.distribution.chi = p.val.distribution.chi,
                       p.val.distribution.chi.scale = p.val.distribution.chi.scale)
  
}

save(results, file = "results_JMI_power.RData")
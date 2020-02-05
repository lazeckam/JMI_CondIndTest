
##### Packages & functions ----

library(bnlearn)

library(Rcpp)
sourceCpp("JMI_components.cpp")

source("functions.R")

##### Grid of parameters size ----

params <- list(gamma = c(1),
               rho = c(0),
               k = 1:10,
               n = c(500, 1000),
               dep = c(TRUE, FALSE))

grid.params <- expand.grid(params)
n.params <- nrow(grid.params)
all.params <- colnames(grid.params)

##### Other parameters size ----

N <- 500
B <- 50

p.val.distribution <- rep(0, N)
p.val.distribution_perm <- rep(0, N)

results <- vector("list", n.params) 

##### Main loop size ----

for(i in 1:n.params){
  
  for(par.name.tmp in all.params) assign(par.name.tmp, grid.params[i, par.name.tmp])
  
  Sigma <- covariance_matrix(k, rho)
  
  cat(paste0(c("|",rep("_", 50), "|", i, "/", n.params), collapse=""), "\n")
  cat("|")
  gwiazdki <- 0
  
  for(j in 1:N){
    
    df <- generate_GTM(k=k, n=n, gamma=gamma, Sigma=Sigma)
    Y <- df$Y
    Z <- as.matrix(df$X[, 2:k])
    
    X <- df$X[, k+1]     # (H_{01}) dep == TRUE  -> X = X_1^{(1)} (M1)
    if(!dep){            # (H_{02}) dep == FALSE -> X = X_1^{(1)} (M2)
      X <- sample(X, n)
    }
    
    cmi_test <- bnlearn::ci.test(as.factor(X), 
                                 as.factor(Y), 
                                 apply(Z, 2, as.factor), test = "mi")
    
    p.val.distribution[j] <- cmi_test$p.value
    
    cmi_perm_test <- bnlearn::ci.test(as.factor(X), 
                                      as.factor(Y), 
                                      apply(Z, 2, as.factor), B=B, test = "sp-mi")
    
    p.val.distribution_perm[j] <- cmi_perm_test$p.value
    
    if(floor(50*j/N - gwiazdki) >= 1){
      replicate(floor(50*j/N - gwiazdki), cat("*"))
      gwiazdki <- gwiazdki+floor(50*j/N - gwiazdki)
    }
    
  }
  cat("|\n")
  
  results[[i]] <- list(parameters = grid.params[i,],
                       p.val.distribution = p.val.distribution,
                       p.val.distribution_perm = p.val.distribution_perm)
  
}

save(results, file = "results_CMI_size.RData")


##### Grid of parameters power ----

params <- list(gamma = c(1),          
               rho = c(0),          
               k = 2:20,            
               n = c(500, 1000),     
               dep = c(TRUE, FALSE)) 
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

##### Other parameters power ----

N <- 500  
B <- 50    

p.val.distribution <- rep(0, N)
p.val.distribution_perm <- rep(0, N)

results <- vector("list", n.params) 

##### Main loop power ----

for(i in 1:n.params){
  
  for(par.name.tmp in all.params) assign(par.name.tmp, grid.params[i, par.name.tmp])
  
  Sigma <- covariance_matrix(k, rho)
  
  cat(paste0(c("|",rep("_", 50), "|", i, "/", n.params), collapse=""), "\n")
  cat("|")
  gwiazdki <- 0
  
  for(j in 1:N){
    
    df <- generate_GTM(k=k, n=n, gamma=gamma, Sigma=Sigma)
    Y <- df$Y
    Z <- as.matrix(df$X[, 2:k])
    
    X <- df$X[, k+1]     # (H_{02}) dep == TRUE  -> X = X_1^{(1)} (M1)
    if(!dep){            # (H_{03}) dep == FALSE -> X = X_1       (M1)
      X <- df$X[, 1] 
    }
    
    cmi_test <- bnlearn::ci.test(as.factor(X), 
                                 as.factor(Y), 
                                 apply(Z, 2, as.factor), test = "mi")
    
    p.val.distribution[j] <- cmi_test$p.value
    
    cmi_perm_test <- bnlearn::ci.test(as.factor(X), 
                                      as.factor(Y), 
                                      apply(Z, 2, as.factor), B=B, test = "sp-mi")
    
    p.val.distribution_perm[j] <- cmi_perm_test$p.value
    
    if(floor(50*j/N - gwiazdki) >= 1){
      replicate(floor(50*j/N - gwiazdki), cat("*"))
      gwiazdki <- gwiazdki+floor(50*j/N - gwiazdki)
    }
    
  }
  cat("|\n")
  
  results[[i]] <- list(parameters = grid.params[i,],
                       p.val.distribution = p.val.distribution,
                       p.val.distribution_perm = p.val.distribution_perm)
  
}
save(results, file = "results_CMI_power.RData")
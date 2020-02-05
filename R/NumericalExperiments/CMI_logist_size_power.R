
##### Packages & functions ----

library(bnlearn)

library(Rcpp)
sourceCpp("JMI_components.cpp")

source("functions.R")

##### Grid of parameters size ----

params <- list(k = c(1:20),
               n = c(500, 1000))

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
  
  p.val.distribution.norm <- rep(NA, N)
  p.val.distribution.chi <- rep(NA, N)
  p.val.distribution.chi.scale <- rep(NA, N)
  
  for(par.name.tmp in all.params) assign(par.name.tmp, grid.params[i, par.name.tmp])
  
  cat(paste0(c("|",rep("_", 50), "|", i, "/", n.params), collapse=""), "\n")
  cat("|")
  gwiazdki <- 0
  
  for(j in 1:N){
    
    df <- generate_logist_0(k=k, n=n)
    Y <- df$Y
    Z <- as.matrix(as.matrix(df$X[, 1:k]))
    X <- df$X[, k+1]
    
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

save(results, file = "results_CMI_size_logist.RData")

##### Grid of parameters power ----

params <- list(k = c(1:20),
               n = c(500, 1000),
               gamma=c(1))

grid.params <- expand.grid(params)
n.params <- nrow(grid.params)
all.params <- colnames(grid.params)

##### Other parameters power ----

N <- 500
B <- 50
p.val.distribution <- rep(0, N)
p.val.distribution_perm <- rep(0, N)

results <- vector("list", n.params) 

##### Main loop power ----

for(i in 1:n.params){
  
  p.val.distribution.norm <- rep(NA, N)
  p.val.distribution.chi <- rep(NA, N)
  p.val.distribution.chi.scale <- rep(NA, N)
  
  for(par.name.tmp in all.params) assign(par.name.tmp, grid.params[i, par.name.tmp])
  
  cat(paste0(c("|",rep("_", 50), "|", i, "/", n.params), collapse=""), "\n")
  cat("|")
  gwiazdki <- 0
  
  for(j in 1:N){
    
    df <- generate_logist_1(k=k, n=n, gamma=gamma)
    Y <- df$Y
    Z <- as.matrix(as.matrix(df$X[, 1:k]))
    X <- df$X[, k+1]
    
    cmi_test <- bnlearn::ci.test(as.factor(X), 
                                 as.factor(Y), 
                                 data.frame(apply(Z, 2, as.factor)), test = "mi")
    
    
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

save(results, file = "results_CMI_power_logist.RData")

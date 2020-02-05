
##### Packages & functions ----

library(Rcpp)
sourceCpp("JMI_components.cpp")

source("functions.R")

##### Grid of parameters ----

params <- list(k = c(1:20),
               n = c(500, 1000))

grid.params <- expand.grid(params)
n.params <- nrow(grid.params)
all.params <- colnames(grid.params)

##### Other parameters ----

N <- 500
B <- 50

choose.distribution <- rep(0, N)
p.val.distribution <- rep(0, N)

results <- vector("list", n.params) 

##### Main loop ----

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
    X <- df$X[, k+1]                   # X = X_1^(1), P(Y=1|X) = ilogit(X_1 + X_2 + ... + X_k), model M3
    
    var.values <- JMI_components(X, Y, Z)
    mean.value <- mean(var.values)
    
    var.mean.full <- replicate(B, {        
      X_perm <- permute(X, Y, Z)$X
      var_values_tmp <- JMI_components(X_perm, Y, Z)
      c(var(var_values_tmp), mean(var_values_tmp))
    })
    
    variance.distribution <- var.mean.full[1,]
    mean.distribution <- var.mean.full[2,]
    var.full <- variance.distribution
    
    var.half <- replicate(B, {             
      X_perm <- permute(X, Y, Z)$X
      i <- sample(1:n, n/2)              
      var(JMI_components(X_perm[i], Y[i], as.matrix(Z[i,])))
    })
    
    test <- t.test(var.full, var.half, "less")
    
    choose.distribution[j] <- test$p.value
    
    if(choose.distribution[j] > 0.05 | mean(2*n*mean.distribution) <= 0){         
      
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

save(results, file = "results_JMI_size_logist.RData")




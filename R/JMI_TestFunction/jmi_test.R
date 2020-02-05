source("jmi_functions.R")
require(bnlearn)

# Test based on JMI statistics
# for testing conditional independence of X and Y given Z

JMI_IndCondTest <- function(X, Y, Z=NULL, alpha=0.05, B=50, alpha.swich=0.05, chi2.dist="chi2.scale"){
  
  # Output: stat:     JMI statistics (value)
  #         rejected: is the null hypothesis rejected (TRUE/FALSE)
  
  stat <- NA
  rejected <- NA
  p.val.distribution <- NA
  
  X <- as.numeric(as.factor(X))
  Y <- as.numeric(as.factor(Y))

  if(is.null(Z)){
    
    test0 <- bnlearn::ci.test(x=as.factor(X) , y=as.factor(Y), test="mi")
    stat <- test0$statistic
    rejected <- (test0$p.value < alpha)
    
  }else{
    
    n <- length(Y)
    Z <- as.matrix(Z)
    
    if (ncol(Z) > 1) {
      Z_one_column <- apply(Z , 1, paste, collapse = "_")
    } else {
      Z_one_column <- Z
    }
    Z_one_column <- factor(Z_one_column)

    var.values <- JMI_components(X, Y, Z)
    
    mean.value <- mean(var.values)
    
    var.mean.full <- matrix(NA, B, 2, dimnames=list(1:B, c("var.dist", "mean.dist")))
    var.half      <- matrix(NA, B, 1, dimnames=list(1:B, c("var.dist")))

    for(b_ind in 1:B){
      
      X_perm <- permute_sample(X, Z_one_column)

      var_values_tmp <- JMI_components(X_perm, Y, Z)
      var.mean.full[b_ind, ] <- c(var(var_values_tmp), mean(var_values_tmp))
      
      X_perm <- permute_sample(X, Z_one_column)
      i <- sample(1:n, n/2)              
      var.half[b_ind, ] <- var(JMI_components(X_perm[i], Y[i], as.matrix(Z[i,])))
      
    }
    
    variance.distribution <- var.mean.full[,1]
    mean.distribution     <- var.mean.full[,2]
    
    var.half <- as.vector(var.half)

    test <- t.test(var.full, var.half, "less")
    
    choose.distribution <- test$p.value
    
    if(choose.distribution > alpha.swich | mean(mean.distribution) <= 0){         
      
      sd.test <- sd(mean.distribution)
      
      p.val.distribution <- pnorm(mean.value, 
                                  mean=mean(mean.distribution), 
                                  sd=sd.test, lower.tail=FALSE)
      
      if(sd.test == 0){
        p.val.distribution <- 1
      }
      
    }else{
      
      if(chi2.dist == "chi2.scale"){
        
        par <-  zhang_params(mean.distribution)
        lambda <- par$d
        alpha  <- par$alpha
        beta   <- par$beta

        p.val.distribution <- pchisq((mean.value - beta) / alpha, 
                                     df=lambda, lower.tail = FALSE)
      }else{
        
        p.val.distribution <- pchisq(mean.value, 
                                     df=mean(mean.distribution), lower.tail = FALSE)
        
      }
    }
    
    stat <- mean.value
    rejected <- (p.val.distribution < alpha)
  }
  
  return(list(stat=stat,
              rejected=rejected))
}

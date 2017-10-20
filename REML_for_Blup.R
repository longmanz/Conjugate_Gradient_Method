#' This is a test script to demonstrat the GREML method in R. 
#' Sept. 21 2017, Longda Jiang
#' 
#' Update: 
#' In this part, I will update the REML method to make it more efficient.
#' version 1.05, on Oct 20, 2017. 




#' 4. REML

log_likelihood <- function(grm, y, sigma_g, sigma_e){
    V <- grm * sigma_g + diag(x = 1, nrow = 500) * sigma_e
    V_i <- solve(V)
    
    # in this example there is no other fixed effect, only the miu (mean), thus:
    X <- rep(1, 500)
    XtVX <- solve(t(X)%*% V_i %*% X)
    
    P <- V_i + V_i %*% X %*% XtVX %*% t(X) %*% V_i
    
    log_value <- -0.5*(log(det(V)) + log(det(XtVX)) + t(y) %*% P %*% y)
    return(log_value)
}


## This is pretty slow...
#  solve2 <- function(x){
#      eigen_list <- eigen(x)
#      eigen_val <- eigen_list$values
#      eigen_vec <- eigen_list$vectors
#      eigen_vec %*% diag(1/eigen_val, nrow = nrow(x)) %*% t(eigen_vec)
#  }
#  
#  system.time(solve(V))
#  system.time(solve2(V))




reml <- function(grm, y, sigma_g=0.95, sigma_e=0.05, iter_num=100, sig_thres = 1e-8){
    

    theta <- rbind(sigma_g, sigma_e)
    X <- rep(1, length(y))
    
    for(i in 1:iter_num){
        
        V <- grm * theta[1] + diag(x = theta[2], nrow = length(y))
        V_i <- solve(V)
        XtViX <- t(X)%*% V_i %*% X
        XtViX_i <- solve(XtViX)
        P <- V_i + V_i %*% X %*% XtViX_i %*% t(X) %*% V_i
        
        ytP <- t(y) %*% P
        Py <- P %*% y
        AP <- grm %*% P
        
        
        ytPAPAy <- ytP %*% AP %*% grm %*% y
        ytPAPPy <- ytP %*% AP %*% Py
        ytPPAPy <- ytP %*% P %*% grm %*% Py
        ytPPPy <- ytP %*% P %*% Py            # they are all simply scalar 
        
        AI <- matrix(c(ytPAPAy, ytPAPPy, ytPPAPy, ytPPPy), nrow=2, byrow = T)
        
        # next is the score function:
        element_1 <- sum(diag(P %*% grm)) - ytP %*% AP %*% y
        element_2 <- sum(diag(P)) - ytP %*% Py
        
        score <- matrix(c(element_1, element_2), nrow=2, byrow = T)
        
        theta2 <- theta - 0.316*solve(AI)%*%score  # this 0.316 here is just adopted from ASREML software
        
        # we need to save the LogL0 for p value calculation
        if(i == 1){
            log_value_0 <- -0.5*(log(det(V)) + log(det(XtViX)) + ytP %*% y)
        }
        
        
        if( abs(theta2[1] -  theta[1]) <= sig_thres){
            theta <- theta2
            log_value <- -0.5*(log(det(V)) + log(det(XtViX)) + ytP %*% y)
            break
        }else{
            theta <- pmax(1e-3, theta2)
            cat(paste("iteration", i, ": sigma_g =", theta[1], "sigma_e =", theta[2], "\n"))
            
            rm(theta2, V, V_i, XtViX_i, P)
            rm(ytP, Py, AP)
            rm( ytPAPAy, ytPAPPy, ytPPAPy, ytPPPy )
            rm(AI, element_1, element_2, score)
        }
    }
    
    # Likelihood ratio test:
    chisq_val <- abs(log_value - log_value_0)*2
    p_val <- pchisq( chisq_val, df=1, lower.tail = F)
    
    list(p_logL=p_val, Vcomp=theta)
    # theta
}

reml(grm, pheno, sig_thres = 1e-5)


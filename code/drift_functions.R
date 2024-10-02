
drift_init <- function(Y, 
                       Kmax = 30, 
                       ebnm_fn = ebnm::ebnm_point_normal, 
                       verbose = FALSE){
  fl <- flash_init(Y) %>%
    flash_set_verbose(1L * verbose) %>%
    flash_greedy(Kmax = Kmax, ebnm_fn = c(ebnm_fn, ebnm::ebnm_normal))
  
  fl_flash_fit <- flash_fit(fl)
  drift_obj <- list(n = nrow(Y),
                    p = ncol(Y),
                    K = fl$n_factors,
                    Y = Y,
                    YYt = tcrossprod(Y),
                    EL = fl$L_pm,
                    EL2 = flash_fit_get_p2m(fl_flash_fit,1),
                    EF = fl$F_pm,
                    EF2 = flash_fit_get_p2m(fl_flash_fit,2),
                    CovF = diag(rep(1, ncol(fl$F_pm))),
                    EFtEF = crossprod(fl$F_pm),
                    resid_s2 = 1/fl$residuals_sd,
                    prior_s2 = sapply(fl$F_ghat, function(x){return(x$sd)})^2,
                    KL_l = fl_flash_fit$KL[[1]],
                    KL_f = sum(fl_flash_fit$KL[[2]]),
                    elbo = fl$elbo,
                    fitted_g_l = fl$L_ghat,
                    fitted_g_F = fl$F_ghat,
                    ebnm_fn = ebnm_fn)
  
  return(drift_obj)
}

# I think the below is correct
normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}

#modified from drift.alpha
# update_elbo <- function(drift_obj) {
#   elbo <- -0.5 * drift_obj$n * drift_obj$p * (1 + log(2*pi)) + sum(drift_obj$KL_l) + drift_obj$KL_f
#   elbo <- elbo - 0.5 * drift_obj$n * drift_obj$p * sum(log(drift_obj$resid_s2))/length(drift_obj$resid_s2)
#   return(elbo)
# }

# this doesn't lead to a change
update_elbo <- function(drift_obj){
#   # assumes constant residual variance
  elbo <- -0.5*drift_obj$n*drift_obj$p*log(2*pi*drift_obj$resid_s2)
  elbo <- elbo - (0.5*(1/drift_obj$resid_s2))*(drift_obj$n * drift_obj$p *(drift_obj$resid_s2)) #this is where the +1 comes from
  elbo <- elbo + sum(drift_obj$KL_l) + drift_obj$KL_f
  return(elbo)
}

drift_update <- function(drift_obj,
                         tol = 1e-6,
                         maxiter = 5000){
  num_iter <- 0
  elbo.list <- c()
  old.elbo <- -Inf
  converged <- FALSE
  print(num_iter)
  
  while ((num_iter < maxiter) & (converged == FALSE)){
  # update posteriors for factors
  ELtL <- crossprod(drift_obj$EL, (drift_obj$EL/drift_obj$resid_s2)) # from drift.alpha code, checked
  diag(ELtL) <- colSums(drift_obj$EL2/drift_obj$resid_s2) # need second moment for diagonal, checked
  drift_obj$CovF <- solve(ELtL + diag(1/drift_obj$prior_s2)) # from drift.alpha code, checked
  
  #drift_obj$EF <- t(t((drift_obj$EL / drift_obj$resid_s2) %*% drift_obj$CovF) %*% drift_obj$Y) # from drift.alpha code, checked
  # doesn't change anything
  drift_obj$EF <- (1/drift_obj$resid_s2)*t((drift_obj$CovF %*% t(drift_obj$EL) %*% drift_obj$Y))
  
  #other stuff that I might need
  prev_EL <- drift_obj$EL 
  
  #below line is copied from drift.alpha, doesn't change results in my tests
  #drift_obj$EFtEF <- drift_obj$CovF %*% t(drift_obj$EL/drift_obj$resid_s2) %*% drift_obj$YYt %*% (drift_obj$EL/drift_obj$resid_s2) %*% drift_obj$CovF 
  drift_obj$EFtEF <- crossprod(drift_obj$EF) # try this
  drift_obj$EF2 <- drift_obj$EF^2 + rep(diag(drift_obj$CovF), each = drift_obj$p) #checked
  
  # update priors for factors
  drift_obj$prior_s2 <- (diag(drift_obj$EFtEF)/drift_obj$p) + diag(drift_obj$CovF) # checked
  
  # update KL for F
  # plug in expression for prior_s2 in KL
  drift_obj$KL_f <- -0.5 * drift_obj$p * (sum(log(drift_obj$prior_s2)) - log(det(drift_obj$CovF))) #checked
  
  # update loadings
  for (k in 1:drift_obj$K){
    s2 <- drift_obj$resid_s2/ (drift_obj$EFtEF[k,k] + drift_obj$p*drift_obj$CovF[k,k]) #checked
    
    # when using these lines from drift.alpha, I get weird results and come across instances where the elbo decreases
    #x_k <- drift_obj$YYt %*% (drift_obj$EL/drift_obj$resid_s2) %*% drift_obj$CovF[,k] # this step confuses me??
    #x_k <- x_k - drift_obj$EL[,-k] %*% drift_obj$EFtEF[-k,k]  # the second part of term (1)
    #x_k <- x_k - drift_obj$p*drift_obj$EL[,-k] %*% drift_obj$CovF[-k,k,drop = FALSE] # term (2) i think
    
    # the below leads to different results than the above
    x_k <- drift_obj$Y %*% drift_obj$EF[, k] - drift_obj$EL[,-k] %*% drift_obj$EFtEF[-k,k] # I think this is right
    x_k <- x_k - drift_obj$p*drift_obj$EL[,-k] %*% drift_obj$CovF[-k,k,drop = FALSE] # term (2) i think
    x <- x_k * s2/drift_obj$resid_s2 # (1/tau_i * sigma^2) factor
    
    e <- drift_obj$ebnm_fn(x = x, s = sqrt(s2), output = ebnm_output_all())
    
    drift_obj$EL[,k] <- e$posterior$mean
    drift_obj$EL2[,k] <- e$posterior$second_moment
    drift_obj$fitted_g_l[[k]] <- e$fitted_g
    drift_obj$KL_l[k] <- e$log_likelihood - normal_means_loglik(x, sqrt(s2), drift_obj$EL[,k], drift_obj$EL2[,k])
  }
  
  # update residual variance parameters
  # (assume constant variance for now) I think this is var_type = 0
  # Need to look at more closely!
  # these lines do not lead to different results in my tests
  #sum1 <- sum(diag(drift_obj$YYt)) - 2 * sum(drift_obj$EL * drift_obj$YYt %*% (prev_EL/drift_obj$resid_s2) %*% drift_obj$CovF)
  #sum1 <- sum1 + sum(crossprod(drift_obj$EL) * drift_obj$EFtEF)
  #sum2 <- sum(colSums(drift_obj$EL2 - drift_obj$EL^2) * (diag(drift_obj$EFtEF) + drift_obj$p*diag(drift_obj$CovF))) # checked
  #sum3 <- drift_obj$p*sum(crossprod(drift_obj$EL) * drift_obj$CovF) # checked
  
  sum1 <- sum((drift_obj$Y - tcrossprod(drift_obj$EL, drift_obj$EF))^2) # this didn't change anything
  sum2 <- sum(colSums(drift_obj$EL2 - drift_obj$EL^2) * (diag(drift_obj$EFtEF) + drift_obj$p*diag(drift_obj$CovF)))
  sum3 <- drift_obj$p*sum(crossprod(drift_obj$EL) * drift_obj$CovF)
  drift_obj$resid_s2 <- (sum1+sum2+sum3)/(drift_obj$n*drift_obj$p) # (1/np) factor
  
  drift_obj$elbo <- update_elbo(drift_obj)
  print(drift_obj$elbo)
  elbo.list[(num_iter)] <- drift_obj$elbo
  converged <- ((drift_obj$elbo - old.elbo) < tol)
  old.elbo <- drift_obj$elbo
  num_iter <- num_iter + 1
  }
  
  return(list(drift_obj = drift_obj, elbo.list = elbo.list, num_iter = num_iter))
}

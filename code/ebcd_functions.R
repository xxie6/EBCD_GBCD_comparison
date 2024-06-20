
normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}

PolarU <- function(A) {
  svdA <- svd(A)
  out <- svdA$u %*% t(svdA$v)
  return(out)
}

ebcd_init <- function(X = NULL,
                      S = NULL,
                      C = NULL,
                      N = NULL) {
  if (!is.null(C)) {
    A <- C
  } else if (!is.null(S)) {
    eigS <- eigen(S)
    C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
    A <- C
  } else if (!is.null(X)) {
    N <- nrow(X)
    
    if (nrow(X) > ncol(X)) {
      S <- Matrix::crossprod(X) / N
      eigS <- eigen(S)
      C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
      A <- C
    } else {
      A <- X
    }
  }
  
  if (is.null(N)) {
    N <- nrow(A)
  }
  
  nrowA <- nrow(A)
  
  tau <- prod(dim(A)) / sum(A^2)
  Z <- matrix(0, nrow = nrowA, ncol = 0)
  EL <- matrix(0, nrow = ncol(A), ncol = 0)
  
  
  ebcd_obj <- list(
    A = A, N = N, nrowA = nrowA,
    tau = tau, Z = Z, EL = EL
  )
  
  
  return(ebcd_obj)
}

#alternative to irlba::irlba is RSpectra::svds()
#irlba is a randomized svd method so I get different results if I don't set a seed
ebcd_greedy <- function(ebcd_obj,
                        Kmax = 1,
                        ebnm_fn = ebnm::ebnm_point_laplace,
                        tol = 1e-6,
                        maxiter = 500) {
  for (K in 1:Kmax) {
    R <- ebcd_obj$A - Matrix::tcrossprod(ebcd_obj$Z, ebcd_obj$EL)
    svd1 <- irlba::irlba(R, nv = 1, nu = 0)
    dv <- svd1$d * svd1$v
    
    l <- dv / sqrt(ebcd_obj$nrowA)
    Rl <- R %*% l
    z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))
    
    ef1 <- list(
      A = R,
      Z = z,
      EL = l,
      maxiter = maxiter,
      ebnm_fn = ebnm_fn,
      N = ebcd_obj$N,
      nrowA = ebcd_obj$nrowA,
      tau = ebcd_obj$tau,
      tol = ebcd_obj$tol
    )
    ef1 <- ebcd_backfit(ef1)
    
    ebcd_obj$EL <- cbind(ebcd_obj$EL, ef1$EL)
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    ebcd_obj$tau <- ef1$tau
  }
  
  ebcd_obj$ebnm_fn <- ebnm_fn
  
  return(ebcd_obj)
}

ebcd_backfit <- function(ebcd_obj,
                         tol = 1e-6,
                         maxiter = 5000) {
  Kmax <- ncol(ebcd_obj$Z)
  ebcd_obj$KL <- rep(0, length = Kmax)
  ebcd_obj$obj.old <- -Inf
  ebcd_obj$vec.obj <- c()
  ebcd_obj$vec.kl <- c() #line I added
  
  for (iter in 1:maxiter) {
    # Shrinkage step
    ebcd_obj$V <- matrix(0, nrow = nrow(ebcd_obj$EL), ncol = ncol(ebcd_obj$EL))
    
    for (k in 1:Kmax) {
      ebnm_fn <- ebcd_obj$ebnm_fn
      x <- c(Matrix::crossprod(ebcd_obj$A, ebcd_obj$Z[, k])) / ebcd_obj$nrowA
      s <- rep(sqrt(1 / (ebcd_obj$N * ebcd_obj$tau)), times=length(x))
      e <- ebnm_fn(x = x, s = s)
      
      ebcd_obj$EL[, k] <- e$posterior$mean
      ebcd_obj$V[, k] <- e$posterior$sd^2
      ebcd_obj$KL[k] <- e$log_likelihood +
        - normal_means_loglik(x, s, ebcd_obj$EL[, k], ebcd_obj$EL[, k]^2 + ebcd_obj$V[, k])
    }
    
    # Rotation step
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    
    # Precision step
    ebcd_obj$tau <- prod(dim(ebcd_obj$A)) / (sum((ebcd_obj$A - ebcd_obj$Z %*% t(ebcd_obj$EL))^2) + ebcd_obj$nrowA * sum(ebcd_obj$V))
    
    # check convergence
    ebcd_obj$obj <- -ebcd_obj$N * ncol(ebcd_obj$A) / 2 * log(2 * pi / ebcd_obj$tau) +
      -(ebcd_obj$N * ebcd_obj$tau / 2) * (
        sum(ebcd_obj$A^2) / ebcd_obj$nrowA - 2 * sum(diag(t(ebcd_obj$A) %*% ebcd_obj$Z %*% t(ebcd_obj$EL))) / ebcd_obj$nrowA + sum(ebcd_obj$EL^2) + sum(ebcd_obj$V)
      ) +
      +sum(ebcd_obj$KL)
    
    ebcd_obj$vec.obj <- c(ebcd_obj$vec.obj, ebcd_obj$obj)
    ebcd_obj$vec.kl <- c(ebcd_obj$vec.kl, sum(ebcd_obj$KL)) #line I added
    #print(ebcd_obj$obj)
    if (iter >= 10 & ebcd_obj$obj - ebcd_obj$obj.old < tol) break
    ebcd_obj$obj.old <- ebcd_obj$obj
  }
  
  return(ebcd_obj)
}

ebcd <- function(X = NULL,
                 S = NULL,
                 C = NULL,
                 N = NULL,
                 Kmax = 5,
                 ebnm_fn = ebnm::ebnm_point_laplace,
                 tol_greedy = 1e-6,
                 maxiter_greedy = 500,
                 tol_backfit = 1e-6,
                 maxiter_backfit = 5000) {
  ebcd_obj <- ebcd_init(X = X, S = S, C = C, N = N) |>
    ebcd_greedy(
      Kmax = Kmax,
      ebnm_fn = ebnm_fn,
      tol = tol_greedy,
      maxiter = maxiter_greedy
    ) |>
    ebcd_backfit(
      tol = tol_backfit,
      maxiter = maxiter_backfit
    )

  return(ebcd_obj)
}

#add onto an already existing initialization to add more factors
#I think the function as is should work
ebcd_greedy_from_object <- function(ebcd_obj,
                        Kmax_add = 1,
                        ebnm_fn = ebnm::ebnm_point_laplace,
                        tol = 1e-6,
                        maxiter = 500) {
  for (K in 1:Kmax_add) {
    R <- ebcd_obj$A - tcrossprod(ebcd_obj$Z, ebcd_obj$EL)
    svd1 <- irlba::irlba(R, nv = 1, nu = 0)
    dv <- svd1$d * svd1$v
    
    l <- dv / sqrt(ebcd_obj$nrowA)
    Rl <- R %*% l
    z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))
    
    ef1 <- list(
      A = R,
      Z = z,
      EL = l,
      maxiter = maxiter,
      ebnm_fn = ebnm_fn,
      N = ebcd_obj$N,
      nrowA = ebcd_obj$nrowA,
      tau = ebcd_obj$tau,
      tol = ebcd_obj$tol
    )
    ef1 <- ebcd_backfit(ef1)
    
    ebcd_obj$EL <- cbind(ebcd_obj$EL, ef1$EL)
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    ebcd_obj$tau <- ef1$tau
  }
  
  #ebcd_obj$ebnm_fn <- ebnm_fn
  
  return(ebcd_obj)
}

ebcd_backfit_from_fit <- function(ebcd_obj,
                         tol = 1e-6,
                         maxiter = 5000) {
  Kmax <- ncol(ebcd_obj$Z)
  #ebcd_obj$KL <- rep(0, length = Kmax)
  #ebcd_obj$obj.old <- -Inf
  #ebcd_obj$vec.obj <- c()
  
  for (iter in 1:maxiter) {
    # Shrinkage step
    ebcd_obj$V <- matrix(0, nrow = nrow(ebcd_obj$EL), ncol = ncol(ebcd_obj$EL))
    
    for (k in 1:Kmax) {
      ebnm_fn <- ebcd_obj$ebnm_fn
      x <- c(crossprod(ebcd_obj$A, ebcd_obj$Z[, k])) / ebcd_obj$nrowA
      s <- rep(sqrt(1 / (ebcd_obj$N * ebcd_obj$tau)), times=length(x))
      e <- ebnm_fn(x = x, s = s)
      
      ebcd_obj$EL[, k] <- e$posterior$mean
      ebcd_obj$V[, k] <- e$posterior$sd^2
      ebcd_obj$KL[k] <- e$log_likelihood +
        - normal_means_loglik(x, s, ebcd_obj$EL[, k], ebcd_obj$EL[, k]^2 + ebcd_obj$V[, k])
    }
    
    # Rotation step
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    
    # Precision step
    ebcd_obj$tau <- prod(dim(ebcd_obj$A)) / (sum((ebcd_obj$A - ebcd_obj$Z %*% t(ebcd_obj$EL))^2) + ebcd_obj$nrowA * sum(ebcd_obj$V))
    
    # check convergence
    ebcd_obj$obj <- -ebcd_obj$N * ncol(ebcd_obj$A) / 2 * log(2 * pi / ebcd_obj$tau) +
      -(ebcd_obj$N * ebcd_obj$tau / 2) * (
        sum(ebcd_obj$A^2) / ebcd_obj$nrowA - 2 * sum(diag(t(ebcd_obj$A) %*% ebcd_obj$Z %*% t(ebcd_obj$EL))) / ebcd_obj$nrowA + sum(ebcd_obj$EL^2) + sum(ebcd_obj$V)
      ) +
      +sum(ebcd_obj$KL)
    
    ebcd_obj$vec.obj <- c(ebcd_obj$vec.obj, ebcd_obj$obj)
    #print(ebcd_obj$obj)
    if (iter >= 10 & ebcd_obj$obj - ebcd_obj$obj.old < tol) break
    ebcd_obj$obj.old <- ebcd_obj$obj
  }
  
  return(ebcd_obj)
}

transform_ebcd_Z <- function(Y, ebcd_obj){
  Y.svd <- svd(Y)
  Y.UV <- Y.svd$u %*% t(Y.svd$v)
  Z_transformed <- Y.UV %*% ebcd_obj$Z
  return(Z_transformed)
}

#initialization for either constant variance or column-wise variance
# ebcd_general_init <- function(X = NULL,
#                       S = NULL,
#                       C = NULL,
#                       N = NULL,
#                       var.type = 0) {
#   if (!is.null(C)) {
#     A <- C
#   } else if (!is.null(S)) {
#     eigS <- eigen(S)
#     C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
#     A <- C
#   } else if (!is.null(X)) {
#     N <- nrow(X)
#     
#     if (nrow(X) > ncol(X)) {
#       S <- crossprod(X) / N
#       eigS <- eigen(S)
#       C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
#       A <- C
#     } else {
#       A <- X
#     }
#   }
#   
#   if (is.null(N)) {
#     N <- nrow(A)
#   }
#   
#   nrowA <- nrow(A)
#   
#   if (var.type == 1){
#     tau <- nrow(A)/colSums(A^2)
#   }
#   else{
#     tau <- prod(dim(A)) / sum(A^2)
#   }
#   Z <- matrix(0, nrow = nrowA, ncol = 0)
#   EL <- matrix(0, nrow = ncol(A), ncol = 0)
#   
#   
#   ebcd_obj <- list(
#     A = A, N = N, nrowA = nrowA,
#     tau = tau, Z = Z, EL = EL
#   )
# }

#alternative to irlba::irlba is RSpectra::svds()
#irlba is a randomized svd method so I get different results if I don't set a seed
# ebcd_colvar_greedy <- function(ebcd_obj,
#                         Kmax = 1,
#                         ebnm_fn = ebnm::ebnm_point_laplace,
#                         tol = 1e-6,
#                         maxiter = 500) {
#   for (K in 1:Kmax) {
#     R <- ebcd_obj$A - tcrossprod(ebcd_obj$Z, ebcd_obj$EL)
#     svd1 <- irlba::irlba(R, nv = 1, nu = 0)
#     dv <- svd1$d * svd1$v
#     
#     l <- dv / sqrt(ebcd_obj$nrowA)
#     Rl <- R %*% l
#     z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))
#     
#     #think about this step
#     ef1 <- list(
#       A = R,
#       Z = z,
#       EL = l,
#       maxiter = maxiter,
#       ebnm_fn = ebnm_fn,
#       N = ebcd_obj$N,
#       nrowA = ebcd_obj$nrowA,
#       tau = ebcd_obj$tau,
#       tol = ebcd_obj$tol
#     )
#     ef1 <- ebcd_colvar_backfit(ef1)
#     
#     ebcd_obj$EL <- cbind(ebcd_obj$EL, ef1$EL)
#     ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
#     ebcd_obj$tau <- ef1$tau
#   }
#   
#   ebcd_obj$ebnm_fn <- ebnm_fn
#   
#   return(ebcd_obj)
# }

#need to check if this runs!
# Note that in practice, one would need to apply some regularization
# ebcd_colvar_backfit <- function(ebcd_obj,
#                          tol = 1e-6,
#                          maxiter = 5000) {
#   Kmax <- ncol(ebcd_obj$Z)
#   ebcd_obj$KL <- rep(0, length = Kmax)
#   ebcd_obj$obj.old <- -Inf
#   ebcd_obj$vec.obj <- c()
#   
#   for (iter in 1:maxiter) {
#     # Shrinkage step
#     ebcd_obj$V <- matrix(0, nrow = nrow(ebcd_obj$EL), ncol = ncol(ebcd_obj$EL))
#     
#     for (k in 1:Kmax) {
#       ebnm_fn <- ebcd_obj$ebnm_fn
#       x <- c(crossprod(ebcd_obj$A, ebcd_obj$Z[, k])) / ebcd_obj$nrowA
#       s <- sqrt(1 / (ebcd_obj$N * ebcd_obj$tau)) #tau is a vector
#       e <- ebnm_fn(x = x, s = s)
#       
#       ebcd_obj$EL[, k] <- e$posterior$mean
#       ebcd_obj$V[, k] <- e$posterior$sd^2
#       ebcd_obj$KL[k] <- e$log_likelihood +
#         - normal_means_loglik(x, s, ebcd_obj$EL[, k], ebcd_obj$EL[, k]^2 + ebcd_obj$V[, k])
#     }
#     
#     # Rotation step
#     ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% diag(ebcd_obj$tau) %*% ebcd_obj$EL)
#     
#     # Precision step
#     ebcd_obj$tau <- ebcd_obj$nrowA / (colSums((ebcd_obj$A - ebcd_obj$Z %*% t(ebcd_obj$EL))^2) + rowSums(ebcd_obj$V))
#     
#     # check convergence
#     ebcd_obj$obj <- -ebcd_obj$N * ncol(ebcd_obj$A) / 2 * log(2 * pi / ebcd_obj$tau) +
#       -(ebcd_obj$N * ebcd_obj$tau / 2) * (
#         sum(ebcd_obj$A^2) / ebcd_obj$nrowA - 2 * sum(diag(t(ebcd_obj$A) %*% ebcd_obj$Z %*% t(ebcd_obj$EL))) / ebcd_obj$nrowA + sum(ebcd_obj$EL^2) + sum(ebcd_obj$V)
#       ) +
#       +sum(ebcd_obj$KL)
#     
#     ebcd_obj$vec.obj <- c(ebcd_obj$vec.obj, ebcd_obj$obj)
#     #print(ebcd_obj$obj)
#     if (iter >= 10 & ebcd_obj$obj - ebcd_obj$obj.old < tol) break
#     ebcd_obj$obj.old <- ebcd_obj$obj
#   }
#   
#   return(ebcd_obj)
# }


# ebcd_colvar <- function(X = NULL,
#                                      S = NULL,
#                                      C = NULL,
#                                      N = NULL,
#                                      Kmax = 5,
#                                      ebnm_fn = ebnm::ebnm_point_laplace,
#                                      tol_greedy = 1e-6,
#                                      maxiter_greedy = 500,
#                                      tol_backfit = 1e-6,
#                                      maxiter_backfit = 5000){
# }







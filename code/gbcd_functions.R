# Initialize the EB-NMF fit to covariance matrix YY' s.t. E[YY'] =
# LL'+ D from an estimate of L without nonnegative constraints.
ax_init_cov_ebnmf <- function(fl, kset = 1:ncol(fl$flash_fit$EF[[1]])) {
  LL <- fl$flash_fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash_fit$EF[[2]][, kset, drop = FALSE]
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash_fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash_fit$EF[[2]][, -kset, drop = FALSE], FF)
  to.keep <- (colSums(LL) > .Machine$double.eps) &
    (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  return(list(LL, FF))
}


# Fit EBMF to covariance matrix YY' s.t. E[YY'] = LL'+ D, where D =
# s2*I and I is an identity matrix.
#
ax_fit_ebmf_to_YY <- function (dat, fl, extrapolate = TRUE, warmstart = TRUE,
                            maxiter = 500, tol = NULL, epsilon = 2e-2,
                            verbose = 1) {
  if (is.matrix(dat) || inherits(dat, "Matrix")) {
    data_diag <- diag(dat)
  } else {
    data_diag <- Matrix::rowSums(dat$U * dat$D * dat$V)
  }
  
  s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
  s2_diff <- Inf
  old_s2 <- 0
  
  ### alternate between estimating s2 and backfitting until convergence
  while(s2 > 0 && abs(s2_diff - 1) > epsilon) {
    if (is.matrix(dat) || inherits(dat, "Matrix")) {
      dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    }
    else{
      dat_minuss2 <- list(U = dat$U, D = dat$D, V = dat$V,
                          S = Matrix::Diagonal(nrow(dat$U), -s2))
    }
    Y2_diff <- sum((data_diag - s2)^2 - (data_diag - old_s2)^2)
    fl <- fl |>
      flash_update_data(dat_minuss2, Y2_diff = Y2_diff) |>
      flash_backfit(extrapolate = extrapolate, warmstart = warmstart,
                    maxiter = maxiter, tol = tol, verbose = verbose)
    
    old_s2 <- s2
    s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat = dat, fl = fl, s2 = s2))
}

# Fit EB-SNMF to Y to estimate F, with L previously calculated from
# GBCD.
#
ax_fit_ebmf_to_Y <- function(Y, fit.cov, corr_thres, maxiter) {
  
  ### scale GEP membership estimates to 0-1 scale, and calculate
  ### Pearson correlations between L and L-tilde
  k.order <- order(fit.cov$pve, decreasing = TRUE)
  fit.L <- fit.cov$L_pm[, k.order]
  fit.L <- t(t(fit.L)/apply(fit.L, 2, max))
  corr <- diag(cor(fit.cov$L_pm[, k.order], fit.cov$F_pm[, k.order]))
  
  ### estimate GEP signatures by fitting EB-SNMF to gene expression
  ### data matrix with fixed L estimated from covariance decomposition
  ### above
  init.F <- t(solve(crossprod(fit.L), as.matrix(Matrix::crossprod(fit.L, Y))))
  fit.snmf <- flash_init(Y, S = 1/sqrt(nrow(Y)), var_type = 2) |>
    flash_factors_init(
      init = list(as.matrix(fit.L), as.matrix(init.F)),
      ebnm_fn = c(ebnm_generalized_binary, ebnm_point_laplace)
    ) |>
    flash_factors_fix(kset = 1:ncol(fit.L), which_dim = "loadings") |>
    flash_backfit(extrapolate = FALSE, maxiter = maxiter, verbose = 1)
  
  ### calculate the z-score and lfsr of GEP signatures by running
  ### linear regression followed by ash
  J <- ncol(Y)
  genes <- colnames(Y)
  F.est <- matrix(0, J, ncol(fit.snmf$L_pm))
  rownames(F.est) <- genes
  colnames(F.est) <- colnames(fit.snmf$L_pm)
  F.se <- F.est
  F.z.ash <- F.est
  F.lfsr.ash <- F.est
  
  for (j in 1:J) {
    y     <- Y[,j]
    dat   <- as.data.frame(cbind(y, fit.snmf$L_pm))
    fit   <- lm(y ~ 0 + ., dat)
    coefs <- summary(fit)$coefficients
    F.est[j,] <- coefs[, "Estimate"]
    F.se[j,]  <- coefs[, "Std. Error"]
  }
  
  for(k in 1:ncol(F.est)){
    fit <- ash(F.est[,k], F.se[,k], mixcompdist = "normal",method = "shrink")
    F.z.ash[,k] <- fit$result$PosteriorMean/fit$result$PosteriorSD
    F.lfsr.ash[,k] <- fit$result$lfsr
  }
  
  ### return the estimated memberships and signatures only for GEPs
  ### whose L and Ltilde from covariance decomposition are highly
  ### concordant
  k.idx <- corr > corr_thres
  L.pm <- fit.snmf$L_pm[, k.idx]
  F.lfc <- fit.snmf$F_pm[, k.idx]/log(2)
  F.z <- F.z.ash[, k.idx]
  F.lfsr <- F.lfsr.ash[, k.idx]
  colnames(L.pm) <- c("Baseline", paste0("GEP", 1:(ncol(L.pm)-1)))
  colnames(F.lfc) <- colnames(L.pm)
  colnames(F.z) <- colnames(L.pm)
  colnames(F.lfsr) <- colnames(L.pm)
  return(list(gbcd_output = list(L = L.pm, F = list(lfc = F.lfc, z_score = F.z, lfsr = F.lfsr)), k.order = k.order, k.idx = k.idx))
}

ax_fit_gbcd <- function (Y, Kmax, prior = ebnm::ebnm_generalized_binary, 
                      maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
                      control = list(), verbose = 1) {
  
  control <- modifyList(ax_fit_gbcd_control_default(), control, keep.null = TRUE)
  extrapolate <- control$extrapolate
  warmstart <- control$warmstart
  corr_thres <- control$corr_thres
  
  ### form the covariance matrix from the cell by gene matrix of gene expression data
  print("Form cell by cell covariance matrix...")
  start_time = proc.time()
  if (2 * ncol(Y) * mean(Y > 0) < nrow(Y)) {
    # Use lowrank plus sparse representation:
    dat <- list(U = Y, D = rep(1 / ncol(Y), ncol(Y)), V = Y)
  } else {
    # Form covariance matrix:
    dat <- Matrix::tcrossprod(Y) / ncol(Y)
  }
  fit.init <- flash_init(dat, var_type = 0)
  runtime = proc.time() - start_time
  print(runtime)
  
  ### fit EBMF with point laplace prior to covariance matrix without considering the diagonal component
  print("Initialize GEP membership matrix L...")
  start_time = proc.time()
  fit.cov <- fit.init |>
    flash_greedy(Kmax = 1, ebnm_fn = ebnm_point_laplace) |>
    flash_greedy(Kmax = Kmax - 1, ebnm_fn = ebnm_point_laplace) |>
    flash_backfit(maxiter = 25, verbose = verbose)
  
  ### fit EBMF with point laplace prior to covariance matrix with the diagonal component
  fit.cov <- ax_fit_ebmf_to_YY(dat = dat, fl = fit.cov, maxiter = maxiter1, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)
  
  ### initialize EB-NMF fit from the EBMF fit with point laplace prior
  print("Estimate GEP membership matrix L...")
  start_time = proc.time()
  cov.init <- ax_init_cov_ebnmf(fit.cov)
  kmax <- which.max(colSums(cov.init[[1]]))
  fit.cov <- fit.init |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, kmax, drop = FALSE]),
      ebnm_fn = ebnm_point_laplace
    ) |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, -c(kmax), drop = FALSE]),
      ebnm_fn = prior
    ) |>
    flash_backfit(extrapolate = FALSE, warmstart = warmstart, maxiter = 25, verbose = verbose)
  
  ### keep at most Kmax factors based on proportion of variance explained and refit EB-NMF to covariance matrix
  kset <- (fit.cov$pve > 0)
  kall <- 1:fit.cov$n_factors
  if(!all(kset)){
    fit.cov <- flash_factors_remove(fit.cov, kset=kall[!kset])
  }
  fit.cov <- ax_fit_ebmf_to_YY(dat = dat, fl = fit.cov, extrapolate = extrapolate, warmstart = warmstart, maxiter = maxiter2, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)
  
  ### compute rescaled version of L
  fit.gbcd.rescale.ldf <- ldf(fit.cov$fl, type = 'i')
  fit.gbcd.rescale.L <- fit.gbcd.rescale.ldf$L %*% diag(sqrt(fit.gbcd.rescale.ldf$D))
  
  ### estimate GEP signatures by fitting EB-SNMF to gene expression data matrix with fixed L estimated from covariance decomposition above
  print("Estimate GEP signature matrix F...")
  start_time = proc.time()
  res0 <- ax_fit_ebmf_to_Y(Y, fit.cov, corr_thres, maxiter3)
  res <- res0$gbcd_output
  runtime = proc.time() - start_time
  print(runtime)
  
  ### only keep columns of rescaled L with high correlation
  fit.gbcd.rescale.L <- fit.gbcd.rescale.L[,res0$k.order][,res0$k.idx]
  
  return(list(res = res, scaled_L = fit.gbcd.rescale.L))
}

ax_fit_gbcd_control_default <- function()
  list(warmstart = TRUE,
       extrapolate = FALSE,
       corr_thres = 0.8)
conditional_mean_sigma_modified <- function(y, x) {
  # When constructing the ll, NA's (i.e. zero sigmas) are not considered in the -sum() 
  # Before if at least one 0 sigma is observed the whole sum  is replaced by NA
  # Additionally, we include nlminb is Nelder-Mead doesn't work
  
  # Starting values and ensure positive fitted standard deviations
  fit1 <- stats::lm(y ~ x - 1)
  fit2 <- stats::lm(abs(fit1$residuals) ~ x - 1)
  fit2$coefficients[1] <- fit2$coefficients[1] - min(0.001, min(fit2$fitted.values))
  b0 <- c(fit1$coefficients, fit2$coefficients)
  
  # Estimate the model under normality
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    -sum(stats::dnorm(x=y, mean=mu, sd=sigma, log=TRUE), na.rm = TRUE)
    #ifelse(all(sigma >= 0), -sum(stats::dnorm(x=y, mean=mu, sd=sigma, log=TRUE)), NA)
  }
  fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="BFGS"), silent=TRUE)
  
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="Nelder-Mead",
                            control=list(maxit=1000000)), silent=TRUE)
  }
  
  if(inherits(fit, "try-error")|| (fit$convergence != 0)) {
    fitforce <- try(suppressWarnings(optimx(b0, function(b) ll(par=b, y=y, x=x), control = list(all.methods=T))), silent=TRUE)  
    fit$par = as.numeric(fitforce["nlminb",1:4])
  }
  b <- fit$par
  
  # Estimated means and standard deviations
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])
  
  list(mu = mu, sigma = sigma)
}


conditional_truncated_variance_modified <- function(y, x, approach) {
  # Call conditional_mean_sigma_modified instead of conditional_mean_sigma
  if (sum(y <= 0) < 2) {
    stop("Not enough negative quantile residuals!")
  }
  
  if (approach == "ind") {
    cv <- rep(stats::var(y[y <= 0]), length(y))
  } else {
    cv <- tryCatch({
      # Get conditional mean and sigma
      mu_sigma <- conditional_mean_sigma_modified(y, x)
      mu <- mu_sigma$mu
      sigma <- mu_sigma$sigma
      
      # Truncated conditional variance
      if (approach == "scl_N") {
        beta <- -mu / sigma
        beta[beta < -30] <- -30
        cv <- sigma^2 * (1 - beta * stats::dnorm(beta)/stats::pnorm(beta) -
                           (stats::dnorm(beta)/stats::pnorm(beta))^2)
      } else if (approach == "scl_sp") {
        # Kernel density estimate of the standardized residuals
        df <- stats::density((y - mu) / sigma, bw = "SJ")
        lower_int_boundary <- min(df$x)
        pdf <- stats::approxfun(df, yleft = 0, yright = 0)
        
        # Truncation points
        beta <- -mu / sigma
        
        # Integrate the truncated pdf by the trapezoidal rule
        div <- 1000
        h <- (max(beta) - lower_int_boundary) / (div - 1)
        b_approx <- seq(lower_int_boundary, max(beta) + h, h)
        midpoint <- b_approx[-div] + h/2
        y0 <- pdf(b_approx)
        y1 <- b_approx * y0
        y2 <- b_approx^2 * y0
        cb <- cumsum(y0[-1] + y0[-div]) / 2 * h
        m1 <- cumsum(y1[-1] + y1[-div]) / 2 * h
        m2 <- cumsum(y2[-1] + y2[-div]) / 2 * h
        cv_approx <- m2 / cb - (m1 / cb)^2
        cv_approx[cv_approx < 0] <- NA
        
        # Approximate the conditional truncated variance
        cv <- sigma^2 * stats::approx(x = midpoint, y = cv_approx, xout = beta)$y
      }
      if (any(is.na(cv)) | any(!is.finite(cv)) | any(cv < 0)) stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit the ", approach, " estimator, switching to the ind approach!"))
      rep(stats::var(y[y <= 0]), length(y))
    })
  }
  
  cv
}

sigma_matrix_modified <- function(object, sigma_est, misspec) {
  # Uses conditional_truncated_variance_modified instead conditional_truncated_variance
  # and cdf_at_quantile_modified instead of cdf_at_quantile
  
  if(!(sigma_est %in% c("ind", "scl_N", "scl_sp")))
    stop("sigma_estimator can be ind, scl_N or scl_sp")
  
  # Extract elements from object
  y <- object$y
  xq <- object$xq
  xe <- object$xe
  coefficients_q <- object$coefficients_q
  coefficients_e <- object$coefficients_e
  alpha <- object$alpha
  n <- length(y)
  kq <- ncol(xq)
  ke <- ncol(xe)
  
  # Transform the data and coefficients
  if (object$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    coefficients_q[1] <- coefficients_q[1] - max_y
    coefficients_e[1] <- coefficients_e[1] - max_y
  }
  
  # Precompute some quantities
  xbq <- as.numeric(xq %*% coefficients_q)
  xbe <- as.numeric(xe %*% coefficients_e)
  uq <- as.numeric(y - xbq)
  
  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xbq, g = "G1_prime", type = object$g1)
  G2_xe <- G_vec(z = xbe, g = "G2", type = object$g2)
  G2_prime_xe <- G_vec(z = xbe, g = "G2_prime", type = object$g2)
  G2_prime_prime_xe <- G_vec(z = xbe, g = "G2_prime_prime", type = object$g2)
  
  # Check the methods in case of sample quantile / es
  if ((kq == 1) & (ke == 1) & sigma_est != "ind") {
    warning("Changed conditional truncated variance estimation to ind!")
    sigma_est <- "ind"
  }
  
  # Estimate the (conditional) truncated variance
  cv <- conditional_truncated_variance_modified(y = uq, x = xq, approach = sigma_est)
  
  # Estimate the CDF at the quantile predictions
  cdf <- cdf_at_quantile_modified(y = y, x = xq, q = xbq)
  
  # Compute sigma
  sigma <- esreg:::sigma_matrix_loop(
    xq = xq, xe = xe, xbq = xbq, xbe = xbe, alpha = alpha,
    G1_prime_xq = G1_prime_xq,
    G2_xe = G2_xe, G2_prime_xe = G2_prime_xe,
    conditional_variance = cv, cdf = cdf,
    include_misspecification_terms = misspec)
  
  sigma
}

vcovA_modified <- function(object, sigma_est = 'scl_N', sparsity = 'nid', misspec = TRUE, bandwidth_estimator = 'Hall-Sheather') {
  # Call lambda_matrix_modified instead lambda_matrix
  # Call sigma_matrix_modified instead sigma_matrix
  lambda <- lambda_matrix_modified(object = object, sparsity = sparsity,
                          bandwidth_estimator = bandwidth_estimator, misspec = misspec)
  lambda_inverse <- solve(lambda)
  sigma <- sigma_matrix_modified(object = object, sigma_est = sigma_est, misspec = misspec)
  n <- length(object$y)
  cov <- 1/n * (lambda_inverse %*% sigma %*% lambda_inverse)
  rownames(cov) <- colnames(cov) <- names(stats::coef(object))
  cov
}

lambda_matrix_modified <- function(object, sparsity, bandwidth_estimator, misspec) {
  # The only difference is that call cdf_at_quantile_modified
  # instead cdf_at_quantile
  if(!(sparsity %in% c("iid", "nid")))
    stop("sparsity can be iid or nid")
  if(!(bandwidth_estimator %in% c("Bofinger", "Chamberlain", "Hall-Sheather")))
    stop("bandwidth_estimator can be Bofinger, Chamberlain or Hall-Sheather")
  
  # Extract elements from object
  y <- object$y
  xq <- object$xq
  xe <- object$xe
  coefficients_q <- object$coefficients_q
  coefficients_e <- object$coefficients_e
  alpha <- object$alpha
  n <- length(y)
  kq <- ncol(xq)
  ke <- ncol(xe)
  
  # Transform the data and coefficients
  if (object$g2 %in% c(1, 2, 3)) {
    max_y <- max(y)
    y <- y - max_y
    coefficients_q[1] <- coefficients_q[1] - max_y
    coefficients_e[1] <- coefficients_e[1] - max_y
  }
  
  # Precompute some quantities
  xbq <- as.numeric(xq %*% coefficients_q)
  xbe <- as.numeric(xe %*% coefficients_e)
  uq <- as.numeric(y - xbq)
  
  # Evaluate G1 / G2 functions
  G1_prime_xq <- G_vec(z = xbq, g = "G1_prime", type = object$g1)
  G1_prime_prime_xq <- G_vec(z = xbq, g = "G1_prime_prime", type = object$g1)
  G2_xe <- G_vec(z = xbe, g = "G2", type = object$g2)
  G2_prime_xe <- G_vec(z = xbe, g = "G2_prime", type = object$g2)
  G2_prime_prime_xe <- G_vec(z = xbe, g = "G2_prime_prime", type = object$g2)
  
  # Check the methods in case of sample quantile / es
  if ((kq == 1) & (ke == 1) & sparsity != "iid") {
    warning("Changed sparsity estimation to iid!")
    sparsity <- "iid"
  }
  
  # Density quantile function
  dens <- density_quantile_function(y = y, x = xq, u = uq, alpha = object$alpha,
                                    sparsity = sparsity, bandwidth_estimator = bandwidth_estimator)
  
  # Conditional CDF evaluated at conditional quantile
  cdf <- cdf_at_quantile_modified(y = y, x = xq, q = xbq)
  
  # Compute lambda
  lambda <- esreg:::lambda_matrix_loop(
    xq = xq, xe = xe, xbq = xbq, xbe = xbe, alpha = alpha,
    G1_prime_xq = G1_prime_xq, G1_prime_prime_xq = G1_prime_prime_xq,
    G2_xe = G2_xe, G2_prime_xe = G2_prime_xe, G2_prime_prime_xe = G2_prime_prime_xe,
    density = dens, cdf = cdf, include_misspecification_terms = misspec)
  
  lambda
}

cdf_at_quantile_modified <- function(y, x, q) {
  # The only difference is that call conditional_mean_sigma_modified(
  # instead of conditional_mean_sigma
  # Get conditional mean and sigma
  mu_sigma <- conditional_mean_sigma_modified(y, x)
  mu <- mu_sigma$mu
  sigma <- mu_sigma$sigma 
  
  # Empirical CDF of standardized data
  cdf <- function(x) stats::ecdf((y - mu) / sigma)(x)
  
  # CDF of standardized quantile predictions
  z <- (q - mu) / sigma
  cdf(z)
}

esr_backtest_modified = function (r, q, e, alpha, version, B = 0, cov_config = list(sparsity = "nid", 
                                                            sigma_est = "scl_N", misspec = TRUE)) {
  
  if (missing(q) & version %in% c(2)) {
    stop("You need to supply VaR forecast `q` for backtest version ", 
         version)
  }
  if (missing(q)) {
    data <- data.frame(r = r, e = e)
  } else {
    data <- data.frame(r = r, q = q, e = e)
  }
  if (version == 1) {
    model <- r ~ e
    h0 <- c(NA, NA, 0, 1)
    one_sided <- FALSE
  } else if (version == 2) {
    model <- r ~ q | e
    h0 <- c(NA, NA, 0, 1)
    one_sided <- FALSE
  } else if (version == 3) {
    model <- I(r - e) ~ e | 1
    h0 <- c(NA, NA, 0)
    one_sided <- TRUE
  } else {
    stop("This is a non-supported backtest version!")
  }
  fit0 <- esreg::esreg(model, data = data, alpha = alpha, g1 = 2, g2 = 1)
  # The original function uses cvovA but conditional_mean_sigma can sometimes do not converge
  # the original function stop() it, but now we use nlminb optimization after BFGS and Nelder-Mead
  # additionally, we force conditional_mean_sigma  to be greater than 0 (ifelse(X>0,<x,0.001))
  cov0 <- vcovA_modified(fit0, sparsity = cov_config$sparsity, 
                       sigma_est = cov_config$sigma_est, misspec = cov_config$misspec,bandwidth_estimator = 'Hall-Sheather')
  s0 <- fit0$coefficients - h0
  mask <- !is.na(h0)
  if (version %in% c(1, 2)) {
    t0 <- as.numeric(s0[mask] %*% solve(cov0[mask, mask]) %*% 
                       s0[mask])
    pv0_1s <- NA
    pv0_2s <- 1 - stats::pchisq(t0, sum(mask))
  }
  else if (version %in% c(3)) {
    t0 <- as.numeric(s0[mask]/sqrt(cov0[mask, mask]))
    pv0_1s <- stats::pnorm(t0)
    pv0_2s <- 2 * (1 - stats::pnorm(abs(t0)))
  }
  if (B > 0) {
    n <- length(r)
    idx <- matrix(sample(1:n, n * B, replace = TRUE), nrow = n)
    bs_estimates <- apply(idx, 2, function(id) {
      tryCatch({
        fitb <- esreg::esreg(model, data = data[id, ], 
                             alpha = alpha, g1 = 2, g2 = 1, early_stopping = 0)
        sb <- fitb$coefficients - fit0$coefficients
        covb <- vcovA_modified(fitb, sparsity = cov_config$sparsity, 
                             sigma_est = cov_config$sigma_est, misspec = cov_config$misspec)
        list(sb = sb, covb = covb)
      }, error = function(e) NA)
    })
    idx_na <- is.na(bs_estimates)
    share_na <- mean(idx_na)
    if (share_na >= 0.05) 
      stop("More than 5% of the bootstrap replications failed!")
    bs_estimates <- bs_estimates[!idx_na]
    if (version %in% c(1, 2)) {
      tb <- sapply(bs_estimates, function(x) {
        as.numeric(x$sb[mask] %*% solve(x$covb[mask, 
                                               mask]) %*% x$sb[mask])
      })
      tb <- tb[!is.na(tb)]
      pvb_2s <- mean(tb >= t0)
      pvb_1s <- NA
    }
    else if (version %in% c(3)) {
      tb <- sapply(bs_estimates, function(x) {
        x$sb[mask]/sqrt(x$covb[mask, mask])
      })
      tb <- tb[!is.na(tb)]
      pvb_2s <- mean(abs(t0) <= abs(tb))
      pvb_1s <- mean(tb <= t0)
    }
  }
  else {
    pvb_2s <- NA
    pvb_1s <- NA
  }
  ret <- list(pvalue_twosided_asymptotic = pv0_2s, pvalue_twosided_bootstrap = pvb_2s)
  if (version %in% c(3)) {
    ret["pvalue_onesided_asymptotic"] <- pv0_1s
    ret["pvalue_onesided_bootstrap"] <- pvb_1s
  }
  ret
}
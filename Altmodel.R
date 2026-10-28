################################################################################
# 7. Alternative Model: Latent Hierarchical Joint Gamma Regression
#
# Same structure as the primary model (Section 2), with gamma margins in place
# of log-normal margins. With Y, Z the latent PSA and PCa progression times and
# X the p-vector of treatment and demographic variables:
#
#   Y = exp(X'beta_Y) * e_Y,   Z = exp(X'beta_Z) * e_Z
#   e_Y ~ Gamma(alpha_Y, alpha_Y),  e_Z ~ Gamma(alpha_Z, alpha_Z)    (mean 1)
#
# so E[Y|X] = exp(X'beta_Y), and Y|X ~ Gamma(shape alpha_Y, rate alpha_Y / exp(X'beta_Y)).
# The PSA/PCa co-pathway dependence enters through a Gaussian copula: with
# U_A = Phi^{-1}(F_A(A)), (U_Y, U_Z) ~ N_2(0, [[1, rho], [rho, 1]]).
# rho plays the same role as rho = sigma_YZ / (sigma_Y sigma_Z) in the primary
# model (the two coincide exactly when the margins are log-normal).
#
# Priors (hyperparameters match Section 3):
#   beta_Y, beta_Z | beta, Sigma_beta  iid N(beta, Sigma_beta),  Sigma_beta = diag(1 / tau_j)
#   tau_j  iid Gamma(a0, b0)       a0 = 2,  b0 = 1
#   beta_j iid N(m0, v0)           m0 = 0,  v0 = 1000
#   alpha_Y, alpha_Z iid Gamma(s0, r0)   s0 = 2, r0 = 1   (mode 1 = exponential; 95% (0.24, 5.6))
#   rho ~ Uniform(-1, 1)
#
# Posterior computation: Metropolis-within-Gibbs, written directly in R because
# JAGS has no way to express a copula between two gamma margins.
#   - latent event times: Gibbs, truncated to each subject's observed interval
#   - (beta_A, log alpha_A), A = Y, Z: adaptive random-walk Metropolis
#   - beta, tau: Gibbs (conjugate, as in the primary model)
#   - rho: random-walk Metropolis on atanh(rho)
#
# gamma_model() mirrors lognorm_model(): same argument order and defaults, one
# chain per call, draws returned as p by 2 by M (betayz) and 2 by 2 by M (Ryz)
# arrays. Ryz is the copula correlation matrix, the counterpart of Sigmayz.
################################################################################

# ==============================================================================
# Data
# ==============================================================================
progression = do(f=\(){
  subdir = "D4320C00015 De-Identified RDB Data"
  vsit = haven::read_sas(paste(getwd(), subdir, "r_visit.sas7bdat",  sep="/"))
  demo = haven::read_sas(paste(getwd(), subdir, "r_dem.sas7bdat",    sep="/"))
  capx = haven::read_sas(paste(getwd(), subdir, "r_caprx.sas7bdat",  sep="/"))
  prog = haven::read_sas(paste(getwd(), subdir, "rd_psa.sas7bdat",   sep="/"))
  secg = haven::read_sas(paste(getwd(), subdir, "r_ecg.sas7bdat",    sep="/"))
  opia = haven::read_sas(paste(getwd(), subdir, "r_opiuse.sas7bdat", sep="/"))
  rcst = haven::read_sas(paste(getwd(), subdir, "rd_rcist.sas7bdat", sep="/"))
  # SUBJ  is subject identifier
  # VISIT is a numerical visit identifier for each subject
  # vsit$VISDYTRT is number of days since treatment
  # demo$ETHGRP   is patient ethnic group
  # demo$AGEGRP   is patient age group
  # capx$CXCLASS  is patient cancer therapy classification
  # capx$BISPHOS  is patient bisphosphonate therapy indicator
  # capx$ANTIANDR is patient anti-androgen therapy indicator
  # prog$PSPCRIT  is indicator of PSA progression
  # secg$ECGEVAL  is indicator of ECG irregularity
  # opia$OPIUSE   is indicator of opiate prescription
  # rcst$RNNEWLES is indicator of new lesions
  # compute full data
  vsit %>% distinct(SUBJ, VISIT, VISDYTRT) %>% # start with patient/visits
    drop_na() %>%
    # join disease progression variables (PSA, ECG, opiate)
    left_join(
      prog %>% distinct(SUBJ, VISIT, PSPCRIT),
      by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
    ) %>%
    left_join(
      secg %>% distinct(SUBJ, VISIT, ECGEVAL),
      by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
    ) %>%
    left_join(
      opia %>% distinct(SUBJ, VISIT, OPIUSE),
      by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
    ) %>%
    left_join(
      rcst %>% distinct(SUBJ, VISIT, RNNEWLES),
      by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
    ) %>%
    # drop visits that did not check for at least one of PSA or ECG or OPIUSE progression
    filter(apply(cbind(PSPCRIT, ECGEVAL, OPIUSE, RNNEWLES), 1, \(r)!all(is.na(r)))) %>%
    arrange(SUBJ, VISIT) %>%
    group_by(SUBJ) %>% # done in preparation for the call to `do()`
    do(\(tib){
      list( # endpoint indicator variable names and the values that indicate progression
        list(name = "PSPCRIT",  prog = "A"),
        list(name = "ECGEVAL",  prog = "1"),
        list(name = "OPIUSE",   prog = "1"),
        list(name = "RNNEWLES", prog = "1")
      ) %>%
        lapply(\(var0){
          tib %>%
            select(SUBJ:VISDYTRT, all_of(var0$name)) %>%
            drop_na(all_of(var0$name)) %>% # drop rows where this outcome isn't measured
            summarise( # compute the outcome interval bounds for all subjects
              l = max(VISDYTRT[VISDYTRT<min(VISDYTRT[.data[[var0$name]]==var0$prog])], 0),
              r = max(min(VISDYTRT[.data[[var0$name]]==var0$prog]), 1),
              .groups="drop"
            ) %>% # update the names in preparation for re-joining
            setColnames(c("SUBJ", paste0(var0$name, c("_lb", "_rb"))))
        }) %>%
        reduce(full_join, by="SUBJ") # join all the outcome intervals together
    }) %>%
    # join demographic variables
    left_join(
      demo %>% distinct(SUBJ, ETHGRP, AGEGRP),
      by=c("SUBJ" = "SUBJ")
    ) %>%
    # join cancer therapy variables
    left_join(
      capx %>% distinct(SUBJ, CXCLASS, BISPHOS, ANTIANDR),
      by=c("SUBJ" = "SUBJ")
    ) %>%
    # clean up demographic and cancer treatment variables
    mutate(
      ETHGRP = as.factor(ETHGRP),
      AGEGRP = as.factor(AGEGRP),
      CXCLASS = as.factor(if_else(CXCLASS=="", "none", as.character(CXCLASS))),
      BISPHOS = as.numeric(BISPHOS),
      ANTIANDR = as.numeric(ANTIANDR)
    )
})
# note: subjects who never progress on an endpoint trigger "no non-missing arguments
# to min; returning Inf" warnings above -- expected, this is what makes _rb = Inf

# design matrix: use the same covariates as the primary log-normal model
x_formula = ~ AGEGRP + ETHGRP + CXCLASS + BISPHOS + ANTIANDR

model_dat = progression %>%
  ungroup() %>%
  filter(!is.na(PSPCRIT_lb), !is.na(RNNEWLES_lb)) %>% # assessed for both outcomes
  drop_na(all_of(all.vars(x_formula)))

x = model.matrix(x_formula, data = model_dat) # n by p, includes intercept

# interval bounds (y = PSA, z = PCa), capped as in Section 4 of the report:
# upper = min(factual upper, 2000), lower = min(factual lower, 1000)
boundsyz = array(
  c(pmin(model_dat$PSPCRIT_lb, 1000), pmin(model_dat$RNNEWLES_lb, 1000),
    pmin(model_dat$PSPCRIT_rb, 2000), pmin(model_dat$RNNEWLES_rb, 2000)),
  dim = c(nrow(x), 2, 2),
  dimnames = list(model_dat$SUBJ, c("y", "z"), c("lower", "upper"))
)

# ==============================================================================
# Model functions
# ==============================================================================
# map times to the copula (standard normal) scale, U = Phi^{-1}(F(t)), stable in both tails
gamma_to_normal = function(t,     # vector of times
                           shape, # gamma shape
                           rate   # gamma rate (vector or scalar)
) {
  lp_lo = pgamma(t, shape, rate, log.p = TRUE)
  lp_hi = pgamma(t, shape, rate, lower.tail = FALSE, log.p = TRUE)
  ifelse(lp_lo < log(0.5), qnorm(lp_lo, log.p = TRUE), -qnorm(lp_hi, log.p = TRUE))
}

# map copula-scale values back to times, t = F^{-1}(Phi(u)), stable in both tails
normal_to_gamma = function(u, shape, rate) {
  ifelse(u < 0,
         qgamma(pnorm(u, log.p = TRUE), shape, rate, log.p = TRUE),
         qgamma(pnorm(u, lower.tail = FALSE, log.p = TRUE), shape, rate,
                lower.tail = FALSE, log.p = TRUE))
}

# draw N(m, s^2) truncated to [a, b], vectorized
rtnorm = function(m, s, a, b) {
  al = (a - m) / s; bl = (b - m) / s
  flip = al > 0 # interval in the upper tail: reflect for precision
  lo = ifelse(flip, -bl, al); hi = ifelse(flip, -al, bl)
  plo = pnorm(lo); phi = pnorm(hi)
  z = qnorm(plo + runif(length(m)) * (phi - plo))
  z = ifelse(is.finite(z), z, ifelse(is.finite(hi), hi, lo))
  z = pmin(pmax(z, lo), hi)
  m + s * ifelse(flip, -z, z)
}

# per-subject gamma log density and copula scores for one outcome
margin_terms = function(t, x, beta_A, alpha_A) {
  rate = alpha_A / exp(drop(x %*% beta_A))
  list(logf = dgamma(t, shape = alpha_A, rate = rate, log = TRUE),
       u    = gamma_to_normal(t, alpha_A, rate))
}

# Gaussian copula log density, summed over subjects
copula_loglik = function(u_y, u_z, rho) {
  r2 = rho^2
  sum(-0.5 * log1p(-r2) - (r2 * (u_y^2 + u_z^2) - 2 * rho * u_y * u_z) / (2 * (1 - r2)))
}

# gamma data-to-posterior model pipeline (one chain)
gamma_model = function(boundsyz,      # n by 2 (y vs z) by 2 (lower vs upper bound) array
                       x,             # n by p design matrix
                       m0 = 0,        # mean of global beta vector entries
                       v0 = 1e+3,     # variance of global beta vector entries
                       a0 = 2,        # shape parameter of betaYZ precision matrix diagonal
                       b0 = 1,        # rate parameter of betaYZ precision matrix diagonal
                       s0 = 2,        # shape parameter of alphaYZ gamma prior
                       r0 = 1,        # rate parameter of alphaYZ gamma prior
                       M = 1e+4,      # number of posterior draws after burn-in
                       burnin = 1e+4, # number of burn-in draws
                       thin = 1,      # post burn-in thinning interval
                       n_adapt = 1e+3 # proposal tuning iterations before burn-in (as jags.model)
) {
  n = nrow(x); p = ncol(x); d = p + 1
  L = boundsyz[, , 1]; U = boundsyz[, , 2]
  exact = (L == U)
  n_iter = n_adapt + burnin + M * thin
  
  # initial values: latent times inside their intervals (random, so chains differ)
  tyz = matrix(L + (U - L) * runif(2 * n, 0.25, 0.75), n, 2)
  tyz = pmax(tyz, 1e-3)
  betayz = matrix(0, p, 2)
  betayz[1, ] = log(colMeans(tyz)) + rnorm(2, 0, 0.2)
  logalpha = rnorm(2, 0, 0.2)
  rho = runif(1, -0.3, 0.3)
  beta = rowMeans(betayz)
  tausqbeta = rep(1, p) # precisions of betayz around beta (as in lognorm_model)
  
  # proposal tuning
  prop_cov = rep(list(diag(0.01, d)), 2)
  hist_theta = rep(list(matrix(NA, n_adapt, d)), 2)
  sd_rho = 0.05; acc_rho_win = 0
  acc = c(0, 0, 0)
  
  # log prior for theta_A = (beta_A, log alpha_A), including the log-alpha Jacobian
  log_prior_theta = function(th) {
    sum(dnorm(th[1:p], beta, 1 / sqrt(tausqbeta), log = TRUE)) +
      dgamma(exp(th[d]), s0, r0, log = TRUE) + th[d]
  }
  
  mt = list(margin_terms(tyz[, 1], x, betayz[, 1], exp(logalpha[1])),
            margin_terms(tyz[, 2], x, betayz[, 2], exp(logalpha[2])))
  
  # storage
  betayz_draws = array(NA, c(p, 2, M))
  alphayz_draws = matrix(NA, 2, M)
  Ryz_draws = array(1, c(2, 2, M))
  keep = 0
  
  for (it in 1:n_iter) {
    adapting = it <= n_adapt
    
    # 1. latent event times, each outcome given the other (Gibbs)
    for (k in 1:2) {
      o = 3 - k
      alpha_k = exp(logalpha[k])
      rate = alpha_k / exp(drop(x %*% betayz[, k]))
      a = gamma_to_normal(L[, k], alpha_k, rate)
      b = gamma_to_normal(U[, k], alpha_k, rate)
      u_new = rtnorm(rho * mt[[o]]$u, sqrt(1 - rho^2), a, b)
      t_new = pmin(pmax(normal_to_gamma(u_new, alpha_k, rate), L[, k]), U[, k])
      tyz[, k] = ifelse(exact[, k], L[, k], pmax(t_new, 1e-6))
      mt[[k]] = margin_terms(tyz[, k], x, betayz[, k], alpha_k)
    }
    
    # 2. (beta_A, log alpha_A), A = Y, Z (adaptive random-walk Metropolis)
    for (k in 1:2) {
      th = c(betayz[, k], logalpha[k])
      th_new = th + drop(t(chol(prop_cov[[k]])) %*% rnorm(d))
      mt_new = margin_terms(tyz[, k], x, th_new[1:p], exp(th_new[d]))
      if (all(is.finite(mt_new$logf)) && all(is.finite(mt_new$u))) {
        cl_old = copula_loglik(mt[[1]]$u, mt[[2]]$u, rho)
        u_y_new = if (k == 1) mt_new$u else mt[[1]]$u
        u_z_new = if (k == 2) mt_new$u else mt[[2]]$u
        cl_new = copula_loglik(u_y_new, u_z_new, rho)
        lr = sum(mt_new$logf) + cl_new + log_prior_theta(th_new) -
          (sum(mt[[k]]$logf) + cl_old + log_prior_theta(th))
        if (log(runif(1)) < lr) {
          th = th_new; mt[[k]] = mt_new
          if (!adapting) acc[k] = acc[k] + 1
        }
      }
      betayz[, k] = th[1:p]; logalpha[k] = th[d]
      if (adapting) {
        hist_theta[[k]][it, ] = th
        if (it >= 200 && it %% 100 == 0) {
          prop_cov[[k]] = (2.38^2 / d) * cov(hist_theta[[k]][(it %/% 2):it, ]) + diag(1e-6, d)
        }
      }
    }
    
    # 3. beta, tausqbeta (Gibbs, conjugate)
    prec = 1 / v0 + 2 * tausqbeta
    beta = rnorm(p, (m0 / v0 + tausqbeta * (betayz[, 1] + betayz[, 2])) / prec, 1 / sqrt(prec))
    tausqbeta = rgamma(p, a0 + 1, b0 + ((betayz[, 1] - beta)^2 + (betayz[, 2] - beta)^2) / 2)
    
    # 4. rho (random-walk Metropolis on atanh(rho), uniform prior on rho)
    rho_new = tanh(atanh(rho) + rnorm(1, 0, sd_rho))
    lr = copula_loglik(mt[[1]]$u, mt[[2]]$u, rho_new) - copula_loglik(mt[[1]]$u, mt[[2]]$u, rho) +
      log1p(-rho_new^2) - log1p(-rho^2)
    accepted = log(runif(1)) < lr
    if (accepted) rho = rho_new
    if (adapting) { # tune toward ~40% acceptance
      acc_rho_win = acc_rho_win + accepted
      if (it %% 100 == 0) {
        sd_rho = sd_rho * ifelse(acc_rho_win / 100 > 0.4, 1.2, 0.85)
        acc_rho_win = 0
      }
    } else if (accepted) acc[3] = acc[3] + 1
    
    # store after adaptation and burn-in, at the thinning interval
    post = it - n_adapt - burnin
    if (post > 0 && post %% thin == 0) {
      keep = keep + 1
      betayz_draws[, , keep] = betayz
      alphayz_draws[, keep] = exp(logalpha)
      Ryz_draws[1, 2, keep] = Ryz_draws[2, 1, keep] = rho
    }
  }
  
  list(
    betayz = betayz_draws,   # p by 2 by M array (counterpart of lognorm_model's betayz)
    alphayz = alphayz_draws, # 2 by M array of gamma shapes
    Ryz = Ryz_draws,         # 2 by 2 by M array (counterpart of lognorm_model's Sigmayz)
    accept_rate = setNames(acc / (burnin + M * thin), c("Y", "Z", "rho")) # aim for 0.15-0.40
  )
}

# sample from gamma model posterior predictive distribution to get M by n by 2 array
gamma_sample_ppd = function(x,             # n by p design matrix
                            betayz_draws,  # p by 2 by M array
                            alphayz_draws, # 2 by M array
                            Ryz_draws      # 2 by 2 by M array
) {
  n = nrow(x)
  M = dim(betayz_draws)[3]
  arrapply(1:M, f=\(m){
    betayz = betayz_draws[,,m]
    alphayz = alphayz_draws[,m]
    L = chol(Ryz_draws[,,m])
    uyz = matrix(rnorm(2*n), ncol=2)%*%L # copula-scale normal scores
    sapply(1:2, \(k){
      normal_to_gamma(uyz[,k], alphayz[k], alphayz[k]/exp(drop(x%*%betayz[,k])))
    })
  })
}

# ==============================================================================
# Fit, diagnostics, and outcome measures
# ==============================================================================

# 4 chains, settings as described in Section 4 of the report
# (set these to whatever the lognormal model's final run uses)
gamma_fits = lapply(1:4, \(ch){
  set.seed(ch)
  gamma_model(boundsyz, x, M = 2500, burnin = 1000, thin = 4)
})
sapply(gamma_fits, \(fit)fit$accept_rate) # Metropolis acceptance rates by chain

# posterior diagnostics: posterior_diagnostics() takes one parameter's draws per
# chain, so run it for each parameter with the chains as separate arguments
p = ncol(x)
gamma_params = c(
  list(rho = \(fit)fit$Ryz[1,2,],
       alphaY = \(fit)fit$alphayz[1,],
       alphaZ = \(fit)fit$alphayz[2,]),
  lapply(1:p, \(j)\(fit)fit$betayz[j,1,]) %>% setNames(paste0("betaY[", colnames(x), "]")),
  lapply(1:p, \(j)\(fit)fit$betayz[j,2,]) %>% setNames(paste0("betaZ[", colnames(x), "]"))
)
gamma_diagnostics = lapply(gamma_params, \(f){
  do.call(posterior_diagnostics, c(lapply(gamma_fits, f), list(lag_max = 50)))
})
sapply(gamma_diagnostics, \(dg)dg$Gelman_Rubin) # want < 1.1
sapply(gamma_diagnostics, \(dg)sum(dg$ESS))     # total effective sample size
gamma_diagnostics$rho$trace_plot
gamma_diagnostics$rho$acf_plot

# pool chains
gamma_draws = list(
  betayz  = array(unlist(lapply(gamma_fits, \(fit)fit$betayz)), c(p, 2, 4 * 2500)),
  alphayz = do.call(cbind, lapply(gamma_fits, \(fit)fit$alphayz)),
  Ryz     = array(unlist(lapply(gamma_fits, \(fit)fit$Ryz)), c(2, 2, 4 * 2500))
)

# posterior summaries (rho quantifies the PSA/PCa co-pathway dependence)
sapply(gamma_params, \(f){
  draws = unlist(lapply(gamma_fits, f))
  c(mean = mean(draws), quantile(draws, c(0.025, 0.975)))
}) %>% t()

# outcome measures at the one-year horizon
ppd_gamma = gamma_sample_ppd(x, gamma_draws$betayz, gamma_draws$alphayz, gamma_draws$Ryz)
cm_gamma  = confusion_matrix(ppd_gamma, t_horizon = 365)
cal_gamma = calib_interval(ppd_gamma, boundsyz)

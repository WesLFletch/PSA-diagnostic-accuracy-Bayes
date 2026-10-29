library("rjags")
library("tidyverse")
library("PipeHelpR") # devtools::install_github("WesLFletch/PipeHelpR")

################################################################################
# 7. Alternative Model: Latent Hierarchical Joint Gamma Regression
#    with a shared latent disease-progression effect
#
# Y, Z are the latent PSA and PCa progression times, X the p-vector of treatment
# and demographic variables, and W_i a subject-level latent disease-progression
# effect shared by both outcomes:
#
#   Y | X, W ~ Gamma(alpha_Y, alpha_Y / exp(X'beta_Y + W))    E[Y | X, W] = exp(X'beta_Y + W)
#   Z | X, W ~ Gamma(alpha_Z, alpha_Z / exp(X'beta_Z + W))    E[Z | X, W] = exp(X'beta_Z + W)
#   W ~ N(0, sigma_W^2),   Y and Z conditionally independent given W
#
# Equivalently log A = X'beta_A + W + log(e_A), e_A ~ Gamma(alpha_A, alpha_A).
# The primary model can be written the same way with normal errors in place of
# log(e_A), so the only change is the subject-level error distribution.
# The implied covariance of (log Y, log Z) given X is
#
#   [ sigma_W^2 + trigamma(alpha_Y)    sigma_W^2                     ]
#   [ sigma_W^2                        sigma_W^2 + trigamma(alpha_Z) ]
#
# returned as Sigmayz, directly comparable with the primary model's Sigma_YZ.
#
# Priors (hyperparameters match Section 3):
#   beta_Y, beta_Z | beta, Sigma_beta  iid N(beta, Sigma_beta),  Sigma_beta = diag(1 / tausqbeta_j)
#   tausqbeta_j iid Gamma(a0, b0)    a0 = 2,  b0 = 1
#   beta_j      iid N(m0, v0)        m0 = 0,  v0 = 1000
#   alpha_Y, alpha_Z iid Gamma(s0, r0)   s0 = 2, r0 = 1   (mode 1 = exponential; 95% (0.24, 5.6))
#   sigma_W ~ Uniform(0, 5)
#
# Fitted with rjags; gamma_model() mirrors lognorm_model() (same arguments and
# return format: betayz p by 2 by M, Sigmayz 2 by 2 by M).
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

# gamma data-to-posterior model pipeline
gamma_model = function(boundsyz,      # n by 2 (y vs z) by 2 (lower vs upper bound) array
                       x,             # n by p design matrix (first column intercept)
                       m0 = 0,        # mean of global beta vector entries
                       v0 = 1e+3,     # variance of global beta vector entries
                       a0 = 2,        # shape parameter of betaYZ precision matrix diagonal
                       b0 = 1,        # rate parameter of betaYZ precision matrix diagonal
                       s0 = 2,        # shape parameter of alphaYZ gamma prior
                       r0 = 1,        # rate parameter of alphaYZ gamma prior
                       M = 1e+4,      # number of posterior draws after burn-in
                       burnin = 1e+4, # number of burn-in draws
                       thin = 1,      # post burn-in thinning interval
                       seed = 1       # JAGS random seed (use a different one per chain)
) {
  model_string = "
  model {
    # data likelihood
    for (i in 1:n) {
      # observations are interval-censored
      for (j in 1:2) { Iyz[i,j] ~ dinterval(yz[i,j], boundsyz[i,j,1:2]) }
      for (j in 1:2) { yz[i,j] ~ dgamma(alphayz[j], rateyz[i,j]) }          # gamma responses
      for (j in 1:2) { rateyz[i,j] = alphayz[j]/exp(muyz[i,j]) }            # so E[yz] = exp(muyz)
      for (j in 1:2) { muyz[i,j] = inprod(x[i,1:p], betayz[1:p,j]) + w[i] } # linear predictor
      w[i] ~ dnorm(0, tauw)                                                 # shared latent progression
    }
    # hierarchical beta coefficients
    for (j in 1:p) { for (k in 1:2) { betayz[j,k] ~ dnorm(beta[j], tausqbeta[j]) } }
    # priors
    beta ~ dmnorm(mubeta[1:p], diagbeta[1:p,1:p])
    for (j in 1:p) { mubeta[j] = m0 }
    for (j in 1:p) { tausqbeta[j] ~ dgamma(a0, b0) } # gamma prior on precision
    for (j in 1:2) { alphayz[j] ~ dgamma(s0, r0) }   # gamma shapes
    sigmaw ~ dunif(0, 5)                             # sd of shared latent effect
    tauw = 1/pow(sigmaw, 2)
  }
  "
  n = nrow(x)
  p = ncol(x)
  datalist = list(
    x = x,
    boundsyz = boundsyz,
    Iyz = matrix(1, nrow=n, ncol=2),
    n = n,
    p = p,
    m0 = m0,
    diagbeta = diag(p)/v0,
    a0 = a0,
    b0 = b0,
    s0 = s0,
    r0 = r0
  )
  yz_init = apply(boundsyz, 1:2, mean) # latent times start inside their intervals
  initslist = list(
    yz = yz_init,
    betayz = rbind(log(colMeans(yz_init)), matrix(0, nrow=p-1, ncol=2)),
    alphayz = c(1, 1),
    sigmaw = 0.5,
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = seed
  )
  out_model = jags.model(textConnection(model_string), data=datalist, inits=initslist)
  update(out_model, n.iter=burnin)
  posterior = coda.samples(
    out_model,
    variable.names=c("betayz", "alphayz", "sigmaw"),
    n.iter=M*thin,
    thin=thin
  )
  betayz = posterior[[1]] %>%
    as_tibble() %>%
    select(starts_with("betayz")) %>%
    mutate(draw = row_number(), .before=everything()) %>%
    pivot_longer(-draw, names_to="var_name", values_to="value") %>%
    mutate(
      j = as.numeric(str_extract_all(var_name, "\\d+") %>% sapply(\(v)v[1])),
      k = as.numeric(str_extract_all(var_name, "\\d+") %>% sapply(\(v)v[2])),
      .before=everything(),
      .keep="unused"
    ) %>%
    xtabs(value~j+k+draw, data=.)
  alphayz = posterior[[1]][, c("alphayz[1]", "alphayz[2]")] %>% as.matrix() %>% t() # 2 by M
  sigmaw = as.vector(posterior[[1]][, "sigmaw"])                                    # length M
  M_out = length(sigmaw)
  # implied covariance of (log Y, log Z) given X, comparable with lognorm_model's Sigmayz
  Sigmayz = sapply(1:M_out, \(m){
    sigmaw[m]^2 + diag(trigamma(alphayz[,m]))
  }) %>%
    array(c(2, 2, M_out))
  list(
    betayz = betayz,   # p by 2 by M array
    alphayz = alphayz, # 2 by M array of gamma shapes
    sigmaw = sigmaw,   # length M vector, sd of shared latent effect
    Sigmayz = Sigmayz  # 2 by 2 by M array
  )
}

# sample from gamma model posterior predictive distribution to get M by n by 2 array
gamma_sample_ppd = function(x,             # n by p design matrix
                            betayz_draws,  # p by 2 by M array
                            alphayz_draws, # 2 by M array
                            sigmaw_draws   # length M vector
) {
  n = nrow(x)
  M = dim(betayz_draws)[3]
  arrapply(1:M, f=\(m){
    betayz = betayz_draws[,,m]
    alphayz = alphayz_draws[,m]
    w = rnorm(n, 0, sigmaw_draws[m]) # shared latent progression for each subject
    muyz = x%*%betayz + w
    sapply(1:2, \(k)rgamma(n, alphayz[k], alphayz[k]/exp(muyz[,k])))
  })
}


# ==============================================================================
# Fit, diagnostics, and outcome measures
# ==============================================================================

# 4 chains, settings as described in Section 4 of the report
# (set these to whatever the lognormal model's final run uses)
gamma_fits = lapply(1:4, \(ch){
  gamma_model(boundsyz, x, M = 2500, burnin = 1000, thin = 4, seed = ch)
})

# correlation of (log Y, log Z), the counterpart of the primary model's rho
rho_draws = \(fit)fit$Sigmayz[1,2,]/sqrt(fit$Sigmayz[1,1,]*fit$Sigmayz[2,2,])

# posterior diagnostics: posterior_diagnostics() takes one parameter's draws per
# chain, so run it for each parameter with the chains as separate arguments
p = ncol(x)
gamma_params = c(
  list(rho = rho_draws,
       sigmaw = \(fit)fit$sigmaw,
       alphaY = \(fit)fit$alphayz[1,],
       alphaZ = \(fit)fit$alphayz[2,]),
  lapply(1:p, \(j)\(fit)as.vector(fit$betayz[j,1,])) %>% setNames(paste0("betaY[", colnames(x), "]")),
  lapply(1:p, \(j)\(fit)as.vector(fit$betayz[j,2,])) %>% setNames(paste0("betaZ[", colnames(x), "]"))
)
gamma_diagnostics = lapply(gamma_params, \(f){
  do.call(posterior_diagnostics, c(lapply(gamma_fits, f), list(lag_max = 50)))
})
sapply(gamma_diagnostics, \(dg)dg$Gelman_Rubin) # want < 1.1
sapply(gamma_diagnostics, \(dg)sum(dg$ESS))     # total effective sample size
gamma_diagnostics$rho$trace_plot
gamma_diagnostics$rho$acf_plot

# posterior summaries
sapply(gamma_params, \(f){
  draws = unlist(lapply(gamma_fits, f))
  c(mean = mean(draws), quantile(draws, c(0.025, 0.975)))
}) %>% t()

# pool chains
n_draws = sum(sapply(gamma_fits, \(fit)length(fit$sigmaw)))
gamma_draws = list(
  betayz  = array(unlist(lapply(gamma_fits, \(fit)as.vector(fit$betayz))), c(p, 2, n_draws)),
  alphayz = do.call(cbind, lapply(gamma_fits, \(fit)fit$alphayz)),
  sigmaw  = unlist(lapply(gamma_fits, \(fit)fit$sigmaw))
)

# outcome measures at the one-year horizon
ppd_gamma = gamma_sample_ppd(x, gamma_draws$betayz, gamma_draws$alphayz, gamma_draws$sigmaw)
cm_gamma  = confusion_matrix(ppd_gamma, t_horizon = 365)
cal_gamma = calib_interval(ppd_gamma, boundsyz)

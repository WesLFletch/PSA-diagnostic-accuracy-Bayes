library("rjags")
library("tidyverse")
library("PipeHelpR") # devtools::install_github("WesLFletch/PipeHelpR")
################################################################################
# 7. Alternative Model: Hierarchical Gamma Model around a Common Latent
#    Progression Time
#
# Y, Z are the latent PSA and PCa progression times and X the p-vector of
# treatment and demographic variables. Each subject has a common latent
# progression time T, the only place covariates enter, and Y, Z scatter around it:
#
#   T | X     ~ Gamma(alpha_T, alpha_T / exp(X'beta))        E[T | X] = exp(X'beta)
#   Y | T     ~ Gamma(alpha_Y, alpha_Y / T)                   E[Y | T] = T
#   Z | T     ~ Gamma(alpha_Z, alpha_Z / (exp(delta) T))      E[Z | T] = exp(delta) T
#   Y and Z conditionally independent given T
#
# delta is the average log-time offset of PCa progression relative to PSA
# progression (Y is centered on T, so only one offset is identifiable).
# Covariate effects are shared: exp(beta_j) is the time ratio for both outcomes.
# This replaces the primary model's hierarchical prior on (beta_Y, beta_Z).
#
# Implied covariance of (log Y, log Z) given X, comparable with Sigma_YZ:
#   [ trigamma(alpha_T) + trigamma(alpha_Y)    trigamma(alpha_T)                     ]
#   [ trigamma(alpha_T)                        trigamma(alpha_T) + trigamma(alpha_Z) ]
#
# Priors:
#   beta_j iid N(m0, v0)                     m0 = 0, v0 = 1000
#   delta ~ N(0, v0)
#   alpha_T, alpha_Y, alpha_Z iid Gamma(s0, r0)   s0 = 2, r0 = 1
#
# Fitted with rjags; gamma_model() follows lognorm_model()'s arguments and style.
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
# gamma common-latent-time data-to-posterior model pipeline
gamma_model = function(boundsyz,      # n by 2 (y vs z) by 2 (lower vs upper bound) array
                       x,             # n by p design matrix (first column intercept)
                       m0 = 0,        # mean of beta vector entries
                       v0 = 1e+3,     # variance of beta vector entries and delta
                       s0 = 2,        # shape parameter of alpha gamma priors
                       r0 = 1,        # rate parameter of alpha gamma priors
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
      # Y and Z scatter around the common latent progression time t
      yz[i,1] ~ dgamma(alphay, alphay/t[i])              # E[Y | t] = t
      yz[i,2] ~ dgamma(alphaz, alphaz/(exp(delta)*t[i])) # E[Z | t] = exp(delta) t
      # common latent progression time, the only place covariates enter
      t[i] ~ dgamma(alphat, alphat/exp(inprod(x[i,1:p], beta[1:p])))
    }
    # priors
    for (j in 1:p) { beta[j] ~ dnorm(m0, 1/v0) }
    delta ~ dnorm(0, 1/v0)
    alphat ~ dgamma(s0, r0)
    alphay ~ dgamma(s0, r0)
    alphaz ~ dgamma(s0, r0)
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
    v0 = v0,
    s0 = s0,
    r0 = r0
  )
  yz_init = apply(boundsyz, 1:2, mean) # latent times start inside their intervals
  initslist = list(
    yz = yz_init,
    t = yz_init[,1],
    beta = c(log(mean(yz_init[,1])), rep(0, p-1)),
    delta = log(mean(yz_init[,2])/mean(yz_init[,1])),
    alphat = 1,
    alphay = 1,
    alphaz = 1,
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = seed
  )
  out_model = jags.model(textConnection(model_string), data=datalist, inits=initslist)
  update(out_model, n.iter=burnin)
  posterior = coda.samples(
    out_model,
    variable.names=c("beta", "delta", "alphat", "alphay", "alphaz"),
    n.iter=M*thin,
    thin=thin
  )
  draws = as.matrix(posterior[[1]])
  alphas = draws[, c("alphat", "alphay", "alphaz")] %>% t() # 3 by M
  M_out = nrow(draws)
  # implied covariance of (log Y, log Z) given X, comparable with lognorm_model's Sigmayz
  Sigmayz = sapply(1:M_out, \(m){
    trigamma(alphas["alphat", m]) + diag(trigamma(alphas[c("alphay", "alphaz"), m]))
  }) %>%
    array(c(2, 2, M_out))
  list(
    beta = draws[, paste0("beta[", 1:p, "]"), drop=FALSE] %>% t(), # p by M array
    delta = draws[, "delta"],                                      # length M vector
    alphas = alphas,                                               # 3 by M (T, Y, Z shapes)
    Sigmayz = Sigmayz                                              # 2 by 2 by M array
  )
}

# sample from gamma model posterior predictive distribution to get M by n by 2 array
gamma_sample_ppd = function(x,           # n by p design matrix
                            beta_draws,  # p by M array
                            delta_draws, # length M vector
                            alpha_draws  # 3 by M array (T, Y, Z shapes)
) {
  n = nrow(x)
  M = ncol(beta_draws)
  arrapply(1:M, f=\(m){
    a = alpha_draws[,m]
    t = rgamma(n, a[1], a[1]/exp(drop(x%*%beta_draws[,m]))) # common latent progression time
    cbind(rgamma(n, a[2], a[2]/t),                          # PSA progression
          rgamma(n, a[3], a[3]/(exp(delta_draws[m])*t)))    # PCa progression
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
       delta = \(fit)fit$delta,
       alphaT = \(fit)fit$alphas[1,],
       alphaY = \(fit)fit$alphas[2,],
       alphaZ = \(fit)fit$alphas[3,]),
  lapply(1:p, \(j)\(fit)fit$beta[j,]) %>% setNames(paste0("beta[", colnames(x), "]"))
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
gamma_draws = list(
  beta   = do.call(cbind, lapply(gamma_fits, \(fit)fit$beta)),
  delta  = unlist(lapply(gamma_fits, \(fit)fit$delta)),
  alphas = do.call(cbind, lapply(gamma_fits, \(fit)fit$alphas))
)

# outcome measures at the one-year horizon
ppd_gamma = gamma_sample_ppd(x, gamma_draws$beta, gamma_draws$delta, gamma_draws$alphas)
cm_gamma  = confusion_matrix(ppd_gamma, t_horizon = 365)
cal_gamma = calib_interval(ppd_gamma, boundsyz)

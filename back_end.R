library("rjags")
library("tidyverse")
library("PipeHelpR") # devtools::install_github("WesLFletch/PipeHelpR")

subdir = "D4320C00015 De-Identified RDB Data"

list.files(paste(getwd(), subdir, sep="/"))

vsit = haven::read_sas(paste(getwd(), subdir, "r_visit.sas7bdat",  sep="/"))
demo = haven::read_sas(paste(getwd(), subdir, "r_dem.sas7bdat",    sep="/"))
prog = haven::read_sas(paste(getwd(), subdir, "rd_psa.sas7bdat",   sep="/"))
secg = haven::read_sas(paste(getwd(), subdir, "r_ecg.sas7bdat",    sep="/"))
opia = haven::read_sas(paste(getwd(), subdir, "r_opiuse.sas7bdat", sep="/"))

# SUBJ  is subject identifier
# VISIT is a numerical visit identifier for each subject
# vsit$VISDYTRT is number of days since treatment
# prog$PSPCRIT  is indicator of PSA progression
# secg$ECGEVAL  is indicator of ECG irregularity
# opia$OPIUSE   is indicator of opiate prescription

# compute full data
progression = vsit %>% distinct(SUBJ, VISIT, VISDYTRT) %>% # start with patient/visits
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
  # drop visits that did not check for at least one of PSA or ECG or OPIUSE progression
  filter(apply(cbind(PSPCRIT, ECGEVAL, OPIUSE), 1, \(r)!all(is.na(r)))) %>%
  arrange(SUBJ, VISIT) %>%
  group_by(SUBJ) %>% # done in preparation for the call to `do()`
  do(\(tib){
    list( # endpoint indicator variable names and the values that indicate progression
      list(name = "PSPCRIT", prog = "A"),
      list(name = "ECGEVAL", prog = "1"),
      list(name = "OPIUSE",  prog = "1")
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
  )

# data-to-posterior pipeline
fit_model = function(boundsyz,       # n by 2 (y vs z) by 2 (lower vs upper bound) array
                     x,              # n by p design matrix
                     m0 = 0,         # mean of global beta vector entries
                     v0 = 1e+3,      # variance of global beta vector entries
                     a0 = 2,         # shape parameter of betaYZ precision matrix diagonal
                     b0 = 1,         # rate parameter of betaYZ precision matrix diagonal
                     d0 = 2,         # df of SigmaYZ Wishart prior
                     W0 = diag(2)/2, # scale matrix of SigmaYZ Wishart prior
                     M = 1e+4,       # number of posterior draws after burn-in
                     burnin = 1e+4,  # number of burn-in draws
                     thin = 1        # post burn-in thinning interval
                     ) {
  model_string = "
  model {
    # data likelihood
    for (i in 1:n) {
      # observations are interval-censored
      for (j in 1:2) { Iyz[i,j] ~ dinterval(yz[i,j], boundsyz[i,j,1:2]) }
      for (j in 1:2) { yz[i,j] = exp(logyz[i,j]) }                    # lognormal responses
      logyz[i,1:2] ~ dmnorm(muyz[i,1:2], Omegayz[1:2,1:2])            # log is norm round lp
      for (j in 1:2) { muyz[i,j] = inprod(x[i,1:p], betayz[1:p,j]) }  # linear predictor
    }
    # hierarchical beta coefficients
    for (j in 1:p) { for (k in 1:2) { betayz[j,k] ~ dnorm(beta[j], tausqbeta[j]) } }
    # priors
    beta ~ dmnorm(mubeta[1:p], diagbeta[1:p,1:p])
    for (j in 1:p) { mubeta[j] = 0 }
    for (j in 1:p) { tausqbeta[j] ~ dgamma(2, 1) }  # equiv to gamma prior on precision
    Omegayz ~ dwish(W0[1:2,1:2], d0)                # Wishart prior on event covariance
    # additional posteriors to be returned
    Sigmayz = inverse(Omegayz)                      # latent outcome covariance matrix
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
    diagbeta = diag(p)/v0,
    W0 = W0,
    d0 = d0
  )
  initslist = list(
    betayz = matrix(rep(0, 2*p), ncol=2),
    Omegayz = diag(2),
    logyz = log(apply(boundsyz, 1:2, mean))
  )
  out_model = jags.model(textConnection(model_string), data=datalist, inits=initslist)
  update(out_model, n.iter=burnin)
  posterior = coda.samples(
    out_model,
    variable.names=c("betayz", "Sigmayz"),
    n.iter=M*thin,
    thin=thin
  )
  list(
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
      xtabs(value~j+k+draw, data=.),
    Sigmayz = posterior[[1]] %>%
      as_tibble() %>%
      select(starts_with("Sigmayz")) %>%
      mutate(draw = row_number(), .before=everything()) %>%
      pivot_longer(-draw, names_to="var_name", values_to="value") %>%
      mutate(
        j = as.numeric(str_extract_all(var_name, "\\d+") %>% sapply(\(v)v[1])),
        k = as.numeric(str_extract_all(var_name, "\\d+") %>% sapply(\(v)v[2])),
        .before=everything(),
        .keep="unused"
      ) %>%
      xtabs(value~j+k+draw, data=.)
  )
}

# compute confusion matrix from posterior and given covariates
confusion_matrix = function(x,              # length p vector, one row from design matrix
                            betayz_draws,   # p by 2 by M array
                            Sigmayz_draws,  # 2 by 2 by M array
                            t_horizon=365   # clinical horizon of interest (days)
                            ) {
  p = length(x)
  M = dim(betayz_draws)[3]
  arrapply(1:M, f=\(m){ # simulate latent Y,Z pair for all posterior draws
    betayz = betayz_draws[,,m]
    Sigmayz = Sigmayz_draws[,,m]
    L = t(chol(Sigmayz))
    muyz = as.vector(x%*%betayz)
    as.vector(exp(muyz+L%*%rnorm(2)))
  }) %>%
    `<`(t_horizon) %>% # compare simulated event times with time threshold
    do(\(mat)2*mat[,1]+mat[,2]) %>% # assign joint outcomes to confusion matrix entries
    table() %>% # tabulate results
    do(\(v){ # ensure tabulation is properly populated
      out = rep(0, 4) %>% setNames(0:3)
      out[names(v)] = v
      out
    }) %>%
    matrix(nrow=2, byrow=T) %>% # reformat into 2x2 matrix
    do(\(tab)tab/sum(tab)) %>% # convert to empirical cell probabilities
    setDimnames(list("Y" = c("0", "1"), "Z" = c("0", "1"))) # rename table entries
}

# check calibration by posterior factual interval inclusion probabilities
calibration = function(boundsyz,     # n by 2 (y vs z) by 2 (lower vs upper bound) array
                       x,            # n by p design matrix
                       betayz_draws, # p by 2 by M array
                       Sigmayz_draws # 2 by 2 by M array
                       ) {
  n = nrow(x)
  p = ncol(x)
  M = dim(betayz_draws)[3]
  arrapply(1:M, f=\(m){ # simulate all obs' Y,Z pair for all posterior draws
    betayz = betayz_draws[,,m]
    Sigmayz = Sigmayz_draws[,,m]
    L = chol(Sigmayz)
    muyz = x%*%betayz
    exp(muyz+matrix(rnorm(2*n), ncol=2)%*%L)
  }) %>% # returns M by n by 2 array
    do(\(arr){ # compute ppd samples' factual interval inclusion indicators
      arrapply(1:M, f=\(m)boundsyz[,,1]<=arr[m,,]&arr[m,,]<=boundsyz[,,2])
    }) %>% # returns M by n by 2 array
    apply(2:3, mean) %>% # take sample inclusion indicator means across ppd samples
    setDimnames(NULL, c("Y", "Z"))
}

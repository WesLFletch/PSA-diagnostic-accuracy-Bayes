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
fit_model = function(boundsyz,
                     x,
                     diagbeta = 1.0e-3*diag(ncol(x)),
                     diagOmegayz = diag(2)/2,
                     dfOmegayz = 2,
                     M = 1e+4,
                     burnin = 1e+4,
                     thin = 1) {
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
    for (j in 1:p) { tausqbeta[j] ~ dgamma(2, 1) }    # equiv to gamma prior on precision
    Omegayz ~ dwish(diagOmegayz[1:2,1:2], dfOmegayz)  # Wishart prior on event covariance
    # additional posteriors to be returned
    Sigmayz = inverse(Omegayz)                        # latent outcome covariance matrix
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
    diagbeta = diagbeta,
    diagOmegayz = diagOmegayz,
    dfOmegayz = dfOmegayz
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
    n.iter=M,
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

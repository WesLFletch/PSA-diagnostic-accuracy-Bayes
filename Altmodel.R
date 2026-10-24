library("tidyverse")
library("PipeHelpR") 

list.files(paste(getwd(), sep="/"))

vsit = haven::read_sas(paste(getwd(), "r_visit.sas7bdat",  sep="/"))
demo = haven::read_sas(paste(getwd(), "r_dem.sas7bdat",    sep="/"))
prog = haven::read_sas(paste(getwd(), "rd_psa.sas7bdat",   sep="/"))
secg = haven::read_sas(paste(getwd(), "r_ecg.sas7bdat",    sep="/"))
opia = haven::read_sas(paste(getwd(), "r_opiuse.sas7bdat", sep="/"))
bonesc = haven::read_sas(paste(getwd(),"r_bonesc.sas7bdat", sep="/"))
rcist = haven::read_sas(paste(getwd(),"rd_rcist.sas7bdat", sep="/"))

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
  left_join(
    bonesc %>% distinct(SUBJ, VISIT, BSMETDIS),
    by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
  ) %>%
  left_join(
    rcist %>% distinct(SUBJ, VISIT, RNNEWLES),
    by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
  ) %>%
  # drop visits that did not check for at least one of PSA or ECG or OPIUSE progression
  filter(apply(cbind(PSPCRIT, ECGEVAL, OPIUSE, BSMETDIS, RNNEWLES), 1, \(r)!all(is.na(r)))) %>%
  arrange(SUBJ, VISIT) %>%
  group_by(SUBJ) %>% # done in preparation for the call to `do()`
  do(\(tib){
    list( # endpoint indicator variable names and the values that indicate progression
      list(name = "PSPCRIT", prog = "A"),
      list(name = "ECGEVAL", prog = "1"),
      list(name = "OPIUSE",  prog = "1"),
      list(name = "BSMETDIS", prog = "1"), 
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
  )

endpoints <- c("PSPCRIT", "ECGEVAL", "OPIUSE", "BSMETDIS", "RNNEWLES")

# 2. Reshape, calculate midpoints, and generate summary statistics
summary_stats_all <- endpoints %>%
  map_dfr(function(endpoint) {
    lb_col <- paste0(endpoint, "_lb")
    rb_col <- paste0(endpoint, "_rb")
    
    progression %>%
      # Select demographics and specific endpoint bounds
      select(SUBJ, ETHGRP, AGEGRP, all_of(c(lb_col, rb_col))) %>%
      # Keep only patients who actually progressed (where right bound is finite/not NA)
      filter(!is.na(.data[[rb_col]]) & is.finite(.data[[rb_col]])) %>%
      mutate(
        Endpoint = endpoint,
        # Calculate the event time midpoint for interval censoring
        Midpoint = (.data[[lb_col]] + .data[[rb_col]]) / 2
      )
  }) %>%
  # 3. Group by Endpoint to calculate summary statistics
  group_by(Endpoint) %>%
  summarize(
    Count        = n(),
    Mean_Midpoint  = mean(Midpoint, na.rm = TRUE),
    Median_Midpoint= median(Midpoint, na.rm = TRUE),
    SD_Midpoint    = sd(Midpoint, na.rm = TRUE),
    Min_Midpoint   = min(Midpoint, na.rm = TRUE),
    Max_Midpoint   = max(Midpoint, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats_all)

library(rjags)
library(coda)
library(tibble)
library(haven)
library(dplyr)
library(tidyr)

# =====================================================================
# 0. HELPERS (same as the lognormal script)
# =====================================================================

# Center/scale varying covariates; add an intercept column if there is none.
standardize_x <- function(x) {
  x <- as.matrix(x)
  varying <- apply(x, 2, sd) > 0
  x[, varying] <- scale(x[, varying, drop = FALSE])
  if (all(varying)) {
    message("No constant column found in x: adding an intercept column.")
    x <- cbind(Intercept = 1, x)
  }
  x
}

# Column names like "betayz[1,2]" laid out column-major
param_cols <- function(name, nr, nc) {
  as.vector(outer(seq_len(nr), seq_len(nc),
                  function(a, b) sprintf("%s[%d,%d]", name, a, b)))
}

# Turn stacked MCMC draws into a d1 x d2 x M array
draws_to_array <- function(draws, name, d1, d2) {
  array(t(draws[, param_cols(name, d1, d2)]), dim = c(d1, d2, nrow(draws)))
}

summarize_draws <- function(d) {
  tibble(
    Measure = names(d),
    Median  = unname(sapply(d, median, na.rm = TRUE)),
    Lower   = unname(sapply(d, quantile, probs = 0.025, na.rm = TRUE)),
    Upper   = unname(sapply(d, quantile, probs = 0.975, na.rm = TRUE))
  )
}

# =====================================================================
# 1. BIVARIATE GAMMA REGRESSION WITH CORRELATED SUBJECT RANDOM EFFECTS
#    Each endpoint: time ~ Gamma(shape_j, rate = shape_j / mean_ij),
#    log(mean_ij) = x_i' beta_j + u_ij, with (u_i1, u_i2) ~ N(0, Sigma_gamma).
#    The random effects are what make the two endpoints dependent.
#    Column 1 / column 2 of boundsyz are the two endpoints (e.g. PSPCRIT, RNNEWLES).
# =====================================================================
gamma_model_string <- "
model {
  for (i in 1:n) {
    u[i, 1:2] ~ dmnorm(mu0[1:2], Omega_u[1:2, 1:2])

    for (j in 1:2) {
      # each endpoint is only known to lie in [lower, upper)
      Iyz[i, j] ~ dinterval(yz[i, j], boundsyz[i, j, 1:2])
      yz[i, j] ~ dgamma(shape[j], rate[i, j])
      rate[i, j] <- shape[j] / mean_yz[i, j]
      log(mean_yz[i, j]) <- inprod(x[i, 1:p], betayz[1:p, j]) + u[i, j]
    }
  }

  # independent normal priors on the coefficients (log-days scale)
  for (j in 1:p) {
    for (k in 1:2) {
      betayz[j, k] ~ dnorm(m0, prec_beta)
    }
  }

  for (j in 1:2) {
    shape[j] ~ dgamma(shape_a, shape_b)
  }

  Omega_u ~ dwish(W0[1:2, 1:2], d0)
  Sigma_gamma <- inverse(Omega_u)
}
"

fit_gamma <- function(boundsyz, x,
                      m0 = 0, v0 = 100,
                      shape_a = 2, shape_b = 0.1,
                      d0 = 3, W0 = diag(2) / 2,
                      n_chains = 3, M = 2000, burnin = 2000, thin = 2) {
  n <- nrow(x)
  p <- ncol(x)
  
  # Floor lower bounds at 0.01 day and cap infinite upper bounds at a large finite value
  lo <- pmax(boundsyz[, , 1], 0.01)
  hi <- pmin(boundsyz[, , 2], 1e7)
  bad <- which(lo >= hi, arr.ind = TRUE)
  if (nrow(bad) > 0) {
    stop("Empty intervals (lower >= upper) at rows: ",
         paste(unique(bad[, 1]), collapse = ", "))
  }
  bounds_use <- array(0, dim = c(n, 2, 2))
  bounds_use[, , 1] <- lo
  bounds_use[, , 2] <- hi
  
  # Starting event times strictly inside every interval
  init_t <- ifelse(hi >= 1e7, pmax(1.5 * lo, lo + 30), (lo + hi) / 2)
  
  data_list <- list(
    x = x, boundsyz = bounds_use,
    Iyz = matrix(1, nrow = n, ncol = 2),   # every observation lies inside its interval
    n = n, p = p,
    mu0 = c(0, 0),
    m0 = m0, prec_beta = 1 / v0,
    shape_a = shape_a, shape_b = shape_b,
    W0 = W0, d0 = d0
  )
  
  make_inits <- function(chain) {
    list(
      betayz  = matrix(rnorm(p * 2, 0, 0.1), nrow = p, ncol = 2),
      shape   = runif(2, 1, 3),
      Omega_u = diag(2),
      yz      = init_t,
      u       = matrix(rnorm(n * 2, 0, 0.1), nrow = n, ncol = 2),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = chain
    )
  }
  
  model <- jags.model(textConnection(gamma_model_string),
                      data = data_list,
                      inits = lapply(seq_len(n_chains), make_inits),
                      n.chains = n_chains)
  update(model, n.iter = burnin)
  coda.samples(model, variable.names = c("betayz", "shape", "Sigma_gamma"),
               n.iter = M * thin, thin = thin)
}

# =====================================================================
# 2. CONVERGENCE DIAGNOSTICS
# =====================================================================
check_convergence <- function(samples, pdf_file = "gamma_trace_plots.pdf") {
  vars <- varnames(samples)
  rhat <- gelman.diag(samples, multivariate = FALSE)$psrf[, "Point est."]
  ess  <- effectiveSize(samples)
  diag_tbl <- tibble(Parameter = vars, Rhat = as.numeric(rhat), ESS = as.numeric(ess))
  
  cat("Worst R-hat values (should be < 1.05):\n")
  print(head(diag_tbl[order(-diag_tbl$Rhat), ], 10))
  cat("\nLowest effective sample sizes (aim for > 400):\n")
  print(head(diag_tbl[order(diag_tbl$ESS), ], 10))
  
  pdf(pdf_file)
  plot(samples)
  dev.off()
  cat("\nTrace plots written to", pdf_file, "\n")
  invisible(diag_tbl)
}

# =====================================================================
# 3. POSTERIOR PREDICTIVE DISTRIBUTION (M x n x 2 array of event times, days)
#    New subjects with the observed covariates: draw random effects from
#    N(0, Sigma_gamma), then gamma noise around the conditional mean.
# =====================================================================
sample_ppd_gamma <- function(x, beta_draws, shape_draws, sigma_draws, idx) {
  n <- nrow(x)
  ppd <- array(0, dim = c(length(idx), n, 2))
  for (r in seq_along(idx)) {
    m  <- idx[r]
    R  <- chol(sigma_draws[, , m])                         # R'R = Sigma_gamma
    u  <- matrix(rnorm(2 * n), ncol = 2) %*% R             # n x 2 random effects
    mu <- exp(x %*% beta_draws[, , m] + u)                 # conditional means
    shp <- matrix(rep(shape_draws[, m], each = n), nrow = n, ncol = 2)
    
    y <- matrix(rgamma(2 * n, shape = as.vector(shp), rate = as.vector(shp / mu)),
                nrow = n, ncol = 2)
    y[is.nan(y)] <- Inf                                    # overflowed means
    ppd[r, , ] <- y
  }
  ppd
}

# =====================================================================
# 4. DIAGNOSTIC ACCURACY AT THE HORIZON, ONE VALUE PER POSTERIOR DRAW
#    test = endpoint used as the diagnostic test, ref = reference endpoint.
# =====================================================================
accuracy_draws <- function(ppd, test, ref, horizon = 365) {
  tpos <- ppd[, , test] < horizon     # M x n
  rpos <- ppd[, , ref]  < horizon
  tibble(
    Sensitivity = rowSums(tpos & rpos)   / rowSums(rpos),
    Specificity = rowSums(!tpos & !rpos) / rowSums(!rpos),
    PPV         = rowSums(tpos & rpos)   / rowSums(tpos),
    NPV         = rowSums(!tpos & !rpos) / rowSums(!tpos),
    Prev_test   = rowMeans(tpos),
    Prev_ref    = rowMeans(rpos)
  )
}

# =====================================================================
# 5. CALIBRATION by interval type
# =====================================================================
calib_by_type <- function(ppd, boundsyz, endpoints) {
  out <- lapply(1:2, function(k) {
    lo <- boundsyz[, k, 1]
    hi <- boundsyz[, k, 2]
    # transpose so the length-n bounds recycle down the n x M matrix correctly
    inside <- (t(ppd[, , k]) >= lo) & (t(ppd[, , k]) <= hi)
    type <- ifelse(is.infinite(hi), "right-censored (upper = Inf)",
                   ifelse(lo == 0 & hi <= 1, "[0, 1]", "finite interval"))
    types <- sort(unique(type))
    tibble(
      Endpoint  = endpoints[k],
      Type      = types,
      N         = as.integer(table(type)[types]),
      Inclusion = as.numeric(tapply(rowMeans(inside), type, mean)[types])
    )
  })
  do.call(rbind, out)
}

# =====================================================================
# 6. CONVENTIONAL (COMPLETE-CASE) ANALYSIS FROM THE OBSERVED BOUNDS
#    progressed by horizon: upper <= horizon;  free at horizon: lower >= horizon;
#    otherwise the status cannot be determined (incomplete follow-up).
# =====================================================================
derive_status <- function(boundsyz, horizon = 365) {
  lo <- boundsyz[, , 1]
  hi <- boundsyz[, , 2]
  status <- matrix(NA_integer_, nrow(lo), ncol(lo))
  status[lo >= horizon] <- 0L
  status[hi <= horizon] <- 1L
  status
}

conventional_accuracy <- function(status, test, ref) {
  cc <- complete.cases(status[, c(test, ref)])
  tp <- sum(status[cc, test] == 1 & status[cc, ref] == 1)
  fn <- sum(status[cc, test] == 0 & status[cc, ref] == 1)
  tn <- sum(status[cc, test] == 0 & status[cc, ref] == 0)
  fp <- sum(status[cc, test] == 1 & status[cc, ref] == 0)
  
  ci <- function(k, n) if (n > 0) binom.test(k, n)$conf.int else c(NA, NA)
  cat("Subjects with both statuses known:", sum(cc), "of", nrow(status), "\n")
  print(matrix(c(tn, fp, fn, tp), 2,
               dimnames = list(Test = c("neg", "pos"), Reference = c("neg", "pos"))))
  
  tibble(
    Measure  = c("Sensitivity", "Specificity", "PPV", "NPV"),
    Estimate = c(tp / (tp + fn), tn / (tn + fp), tp / (tp + fp), tn / (tn + fn)),
    Lower    = c(ci(tp, tp + fn)[1], ci(tn, tn + fp)[1], ci(tp, tp + fp)[1], ci(tn, tn + fn)[1]),
    Upper    = c(ci(tp, tp + fn)[2], ci(tn, tn + fp)[2], ci(tp, tp + fp)[2], ci(tn, tn + fn)[2])
  )
}

# =====================================================================
# 7. DATA PREPARATION FROM THE TRIAL FILES
#    rd_psa has one row per subject per visit, but PSPCRIT, PSPDYRND (day of
#    progression or censoring, from randomization) and PSPCENS describe the
#    SUBJECT's PSA time-to-event outcome and repeat on every row. So PSA
#    progression is taken from those subject-level variables, not from
#    visit-by-visit PSPCRIT values.
#      PSPCRIT: "A" = PSA progression, "0" = no progression,
#               "5" = death without progression (treated as censored at that day)
#    RNNEWLES (rd_rcist, one value per visit): "1" = new lesions
# =====================================================================
load_trial_visits <- function() {
  # the .sas7bdat files are read straight from the working directory
  rd <- function(f) haven::read_sas(f)
  vsit <- rd("r_visit.sas7bdat")      # SUBJ, VISIT, VISDYTRT (days since treatment)
  prog <- rd("rd_psa.sas7bdat")       # PSPCRIT, PSPDYRND, PSPCENS
  rcst <- rd("rd_rcist.sas7bdat")     # RNNEWLES
  demo <- rd("r_dem.sas7bdat")        # ETHGRP, AGEGRP
  
  visits <- vsit %>%
    distinct(SUBJ, VISIT, VISDYTRT) %>%
    drop_na() %>%
    left_join(distinct(rcst, SUBJ, VISIT, RNNEWLES), by = c("SUBJ", "VISIT")) %>%
    mutate(VISDYTRT = as.numeric(VISDYTRT),
           RNNEWLES = as.character(RNNEWLES)) %>%
    filter(!is.na(RNNEWLES)) %>%                     # keep visits with a lesion assessment
    arrange(SUBJ, VISIT)
  
  dups <- visits %>% count(SUBJ, VISIT) %>% filter(n > 1)
  if (nrow(dups) > 0) {
    warning(nrow(dups), " SUBJ/VISIT pairs appear more than once (conflicting RNNEWLES records).")
  }
  list(visits = visits, prog = prog, demo = distinct(demo, SUBJ, ETHGRP, AGEGRP))
}

# One row per subject: PSA outcome code and the day it happened / was censored
psa_subject_times <- function(prog) {
  prog %>%
    mutate(PSPCRIT = as.character(PSPCRIT)) %>%
    group_by(SUBJ) %>%
    summarise(
      n_crit = n_distinct(PSPCRIT),
      n_day  = n_distinct(PSPDYRND[!is.na(PSPDYRND)]),
      crit   = if (any(PSPCRIT == "A", na.rm = TRUE)) "A" else first(na.omit(PSPCRIT)),
      day    = if (all(is.na(PSPDYRND))) NA_real_ else min(PSPDYRND, na.rm = TRUE),
      .groups = "drop"
    )
}

# Checks that PSPCRIT / PSPDYRND really are constant within a subject
diagnose_psa <- function(prog, psa) {
  cat("PSPCRIT by PSPCENS (all rd_psa rows):\n")
  print(table(PSPCRIT = prog$PSPCRIT, PSPCENS = prog$PSPCENS, useNA = "ifany"))
  cat("\nSubjects in rd_psa:", nrow(psa), "\n")
  cat("Subjects with more than one PSPCRIT value across their rows:", sum(psa$n_crit > 1), "\n")
  cat("Subjects with more than one PSPDYRND value across their rows:", sum(psa$n_day > 1), "\n")
  cat("Subjects with no PSPDYRND:", sum(is.na(psa$day)), "\n")
  cat("\nPSPCRIT per subject:\n")
  print(table(psa$crit, useNA = "ifany"))
  cat("\nPSPDYRND (days from randomization) by PSPCRIT:\n")
  print(tapply(psa$day, psa$crit, summary))
  cat("\n")
  invisible(psa)
}

# Interval bounds for the RNNEWLES endpoint, one row per subject:
#   event seen   -> (last visit day before the first event, first event day)
#   no event     -> (last visit day, Inf)   i.e. right-censored
# Lower bound is at least 0 and upper bound is at least 1.
subject_bounds <- function(visits, outcome, event_value) {
  visits %>%
    filter(!is.na(.data[[outcome]])) %>%
    group_by(SUBJ) %>%
    summarise(
      event_day = if (any(.data[[outcome]] == event_value)) {
        min(VISDYTRT[.data[[outcome]] == event_value])
      } else Inf,
      lb = max(VISDYTRT[VISDYTRT < event_day], 0),
      ub = max(event_day, 1),
      .groups = "drop"
    ) %>%
    select(SUBJ, lb, ub)
}

# Build the bounds array (n x 2 x 2), design matrix and endpoint labels.
make_analysis_data <- function() {
  tr  <- load_trial_visits()
  psa <- psa_subject_times(tr$prog)
  diagnose_psa(tr$prog, psa)
  
  # PSPCRIT: progression on day t -> one-day interval (t - 1, t];
  #          no progression (or death without progression) -> censored at day t
  b_psa <- psa %>%
    filter(!is.na(day), !is.na(crit)) %>%
    mutate(event = crit == "A",
           PSPCRIT_lb = ifelse(event, pmax(day - 1, 0), pmax(day, 0)),
           PSPCRIT_ub = ifelse(event, pmax(day, 1), Inf)) %>%
    select(SUBJ, PSPCRIT_lb, PSPCRIT_ub)
  
  b_rad <- subject_bounds(tr$visits, "RNNEWLES", "1") %>%
    rename(RNNEWLES_lb = lb, RNNEWLES_ub = ub)
  
  cat("Subjects with a PSA outcome:", nrow(b_psa),
      "| with lesion assessments:", nrow(b_rad), "\n")
  
  dat <- b_psa %>%
    inner_join(b_rad, by = "SUBJ") %>%
    inner_join(tr$demo, by = "SUBJ")
  n_before <- nrow(dat)
  dat <- drop_na(dat)
  cat("Subjects with both endpoints and covariates:", nrow(dat),
      "(", n_before - nrow(dat), "dropped for missing covariates )\n")
  
  x <- model.matrix(~ factor(ETHGRP) + factor(AGEGRP), data = dat)
  
  bounds <- array(0, dim = c(nrow(dat), 2, 2))
  bounds[, 1, 1] <- dat$PSPCRIT_lb
  bounds[, 1, 2] <- dat$PSPCRIT_ub
  bounds[, 2, 1] <- dat$RNNEWLES_lb
  bounds[, 2, 2] <- dat$RNNEWLES_ub
  
  list(bounds = bounds, x = x, subj = dat$SUBJ,
       endpoints = c("PSPCRIT", "RNNEWLES"))
}

# =====================================================================
# RUN
# =====================================================================
analysis <- make_analysis_data()

boundsyz_data <- analysis$bounds
x_data        <- analysis$x
endpoints     <- analysis$endpoints

horizon  <- 365    # 12 months in days
test_idx <- 1      # PSA progression (PSPCRIT) is the test, per the project prompt
ref_idx  <- 2      # radiographic progression (RNNEWLES) is the reference
# To flip the direction (RNNEWLES predicting PSPCRIT), swap test_idx and ref_idx.

x_std <- standardize_x(x_data)

samples <- fit_gamma(boundsyz_data, x_std)
check_convergence(samples)

draws       <- as.matrix(samples)                       # all chains stacked
p           <- ncol(x_std)
beta_draws  <- draws_to_array(draws, "betayz",      p, 2)
sigma_draws <- draws_to_array(draws, "Sigma_gamma", 2, 2)
shape_draws <- t(draws[, c("shape[1]", "shape[2]")])    # 2 x M

set.seed(1)
idx <- sample(nrow(draws), min(2000, nrow(draws)))
ppd <- sample_ppd_gamma(x_std, beta_draws, shape_draws, sigma_draws, idx)

# Model-based accuracy with posterior uncertainty
acc <- accuracy_draws(ppd, test_idx, ref_idx, horizon)
cat("\nModel-based accuracy at", horizon, "days (posterior median, 95% CrI):\n")
print(summarize_draws(acc))

# Does the model reproduce the observed data?
cat("\nCalibration by interval type:\n")
print(calib_by_type(ppd, boundsyz_data, endpoints))

status <- derive_status(boundsyz_data, horizon)
cat("\nObserved 12-month status (1 = progressed, 0 = free, NA = undetermined):\n")
for (k in 1:2) {
  cat(endpoints[k], ": ", sum(status[, k] == 1, na.rm = TRUE), " progressed, ",
      sum(status[, k] == 0, na.rm = TRUE), " free, ",
      sum(is.na(status[, k])), " undetermined\n", sep = "")
}

# Conventional complete-case analysis for comparison
cat("\nConventional complete-case analysis:\n")
conv <- conventional_accuracy(status, test_idx, ref_idx)
print(conv)
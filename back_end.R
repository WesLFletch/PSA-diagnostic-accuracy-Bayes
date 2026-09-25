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

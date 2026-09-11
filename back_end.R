setwd("C:/Users/WesLF/Desktop/TAMU/Semesters/2026 Fall/STAT-638-600/Project")
library("tidyverse")
library("PipeHelpR") # devtools::install_github("WesLFletch/PipeHelpR")

subdir = "D4320C00015 De-Identified RDB Data"

list.files(paste(getwd(), subdir, sep="/"))

vsit = haven::read_sas(paste(getwd(), subdir, "r_visit.sas7bdat", sep="/"))
prog = haven::read_sas(paste(getwd(), subdir, "rd_psa.sas7bdat", sep="/"))
cist = haven::read_sas(paste(getwd(), subdir, "rd_rcist.sas7bdat", sep="/"))

# SUBJ is subject identifier
# VISIT is a numerical visit identifier for each subject
# VISDYTRT is number of days since treatment
# prog$PSPCRIT is indicator of overall progression
# cist$OVBRESP is indicator of radiographic metastatic progression

# a primitive join attempt
psa = vsit %>%
  distinct(SUBJ, VISIT, VISDYTRT) %>%
  full_join(
    prog %>% distinct(SUBJ, VISIT, PSPCRIT),
    by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
  ) %>%
  full_join(
    cist %>% distinct(SUBJ, VISIT, OVBRESP),
    by=c("SUBJ" = "SUBJ", "VISIT" = "VISIT")
  ) %>%
  drop_na(SUBJ, VISIT, VISDYTRT) %>%
  filter(apply(cbind(PSPCRIT, OVBRESP), 1, \(r)!all(is.na(r)))) %>%
  replace_na(list("PSPCRIT" = "0", "OVBRESP" = "6")) %>%
  arrange(SUBJ, VISIT)

# get interval bounds of overall progression and radiographic metastatic progression
psa_intervals = psa %>%
  group_by(SUBJ) %>%
  summarise(
    op_ub = min(VISDYTRT[PSPCRIT!="0"]),
    op_lb = max(VISDYTRT[VISDYTRT<min(VISDYTRT[PSPCRIT!="0"])]),
    rmp_ub = min(VISDYTRT[OVBRESP=="3"]),
    rmp_lb = max(VISDYTRT[VISDYTRT<min(VISDYTRT[OVBRESP=="3"])]),
    .groups="drop"
  )

# how do the bounds relate to the 1-year horizon?
mean(psa_intervals$op_ub<365)
mean(psa_intervals$op_lb>365)
mean(psa_intervals$rmp_ub<365)
mean(psa_intervals$rmp_lb>365)

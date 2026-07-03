# descr_LSAC
#
# This script will clean and describe the LSAC data used for the analysis.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
LSAC <- data.frame(read_dta("../data/LSAC.dta"))

# Clean data 
LSAC <- dplyr::mutate(LSAC, 
                      sex = haven::as_factor(sex),
                      indstat = haven::as_factor(indstat),
                      noneng = haven::as_factor(noneng),
                      overwt_w1 = haven::as_factor(overwt_w1),
                      poor_hrqol_w4= haven::as_factor(low_hrqol_w4)) 

# Add missingness indicator for HRQoL
LSAC$m_hrqol_w4 <- ifelse(is.na(LSAC$hrqol_w4), 1, 0)

# Set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# Describe missingness in all analysis variables
# Variable summaries overall
LSAC[, c(ana_vars, "hrqol_w4", "low_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         missing_text = "missing") 

# Reduce analysis sample so that exposure and covariates are complete
LSAC <- LSAC[-which(rowSums(1*is.na(LSAC[, c("sex", "indstat", "noneng", "sep_w1", "overwt_w1")])) >= 1), ]

# Convert age from months to years
LSAC$age_w1 <- LSAC$age_w1 / 12

# Variable summaries overall
LSAC[, c(ana_vars, "hrqol_w4", "low_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                         all_dichotomous() ~ "{n} ({p})"),
                         missing_text = "missing") 

# Variable summaries by exposure
LSAC[, c(ana_vars, "hrqol_w4", "poor_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                         all_dichotomous() ~ "{n} ({p})"),
                         by = "overwt_w1",
                         missing_text = "missing") 

# Missing data patterns
miss_pattern <- mice::md.pattern(LSAC[, c(ana_vars, "hrqol_w4", "poor_hrqol_w4")], 
                                 rotate.names = TRUE)
# Complete cases = 4056

# Variable summaries by missingness in HRQoL
LSAC[, c(ana_vars, "hrqol_w4", "poor_hrqol_w4", "m_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         by = "m_hrqol_w4",
                         missing_text = "missing")

# Save clean data
save(LSAC, file = "../data/LSAC_AD.dta")

# Remove objects from global environment
rm("LSAC", "miss_pattern", "ana_vars")

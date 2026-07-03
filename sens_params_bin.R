# sens_params_bin
#
# This script will estimate relationships between poor HRQoL at age 10-11 years 
# and earlier measures of HRQoL. It is used to inform the choice of sensitivity 
# parameter values for the analysis of the binary outcome.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
load(file = "../data/LSAC_AD.dta")

# Define set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# Create poor HRQoL variables for waves 1 to 3
LSAC$poor_hrqol_w1 <- factor(ifelse(LSAC$hrqol_w1 < 50, "HRQoL < 50", "HRQoL >= 50"),
                            levels = c("HRQoL >= 50", "HRQoL < 50"))
LSAC$poor_hrqol_w2 <- factor(ifelse(LSAC$hrqol_w2 < 50, "HRQoL < 50", "HRQoL >= 50"),
                            levels = c("HRQoL >= 50", "HRQoL < 50"))
LSAC$poor_hrqol_w3 <- factor(ifelse(LSAC$hrqol_w3 < 50, "HRQoL < 50", "HRQoL >= 50"),
                            levels = c("HRQoL >= 50", "HRQoL < 50"))

# Extract relevant variables
LSAC <- LSAC[, c("poor_hrqol_w1", "poor_hrqol_w2", "poor_hrqol_w3", "poor_hrqol_w4", 
                 ana_vars)]

# Add missingness indicators
LSAC$m_poor_hrqol_w1 <- ifelse(is.na(LSAC$poor_hrqol_w1), 1, 0)
LSAC$m_poor_hrqol_w2 <- ifelse(is.na(LSAC$poor_hrqol_w2), 1, 0)
LSAC$m_poor_hrqol_w3 <- ifelse(is.na(LSAC$poor_hrqol_w3), 1, 0)
LSAC$m_poor_hrqol_w4 <- ifelse(is.na(LSAC$poor_hrqol_w4), 1, 0)

# Examine missing data patterns
miss_pattern <- mice::md.pattern(LSAC[, c("poor_hrqol_w1", "poor_hrqol_w2", 
                                          "poor_hrqol_w3", "poor_hrqol_w4")], 
                                 rotate.names = TRUE)

# Poor HRQoL summaries overall
LSAC[, c("poor_hrqol_w1", "m_poor_hrqol_w1", "poor_hrqol_w2", "m_poor_hrqol_w2", 
         "poor_hrqol_w3", "m_poor_hrqol_w3")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         missing_text = "missing")

# Poor HRQoL summaries by missingness in poor HRQoL at wave 4
LSAC[, c("poor_hrqol_w1", "m_poor_hrqol_w1", "poor_hrqol_w2", "m_poor_hrqol_w2", 
         "poor_hrqol_w3", "m_poor_hrqol_w3", "poor_hrqol_w4", "m_poor_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         by = "m_poor_hrqol_w4",
                         missing_text = "missing") 

# Fit univariable logistic regression models to investigate odds ratio for poor 
# hrqol at earlier waves between those with and without missing hrqol at wave 4
logit_w1 <- glm(poor_hrqol_w1 ~ m_poor_hrqol_w4, family = "binomial", data = LSAC)
logit_w2 <- glm(poor_hrqol_w2 ~ m_poor_hrqol_w4, family = "binomial", data = LSAC)
logit_w3 <- glm(poor_hrqol_w3 ~ m_poor_hrqol_w4, family = "binomial", data = LSAC)
exp(coef(logit_w1)[2])
exp(coef(logit_w2)[2])
exp(coef(logit_w3)[2])

# Fit multivariable logistic regression models to investigate log odds ratio for 
# poor hrqol at earlier waves between those with and without missing hrqol at wave 4
glm(poor_hrqol_w1 ~ m_poor_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)
glm(poor_hrqol_w2 ~ m_poor_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)
glm(poor_hrqol_w3 ~ m_poor_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)

# Fit univariable logistic regression models to investigate odds ratio for 
# missingness in hrqol at wave 4 between those with and without poor hrqol at 
# earlier waves - sane as above
logit_2_w1 <- glm(m_poor_hrqol_w4 ~ poor_hrqol_w1, family = "binomial", data = LSAC)
logit_2_w2 <- glm(m_poor_hrqol_w4 ~ poor_hrqol_w2, family = "binomial", data = LSAC)
logit_2_w3 <- glm(m_poor_hrqol_w4 ~ poor_hrqol_w3, family = "binomial", data = LSAC)
exp(coef(logit_2_w1)[2])
exp(coef(logit_2_w2)[2])
exp(coef(logit_2_w3)[2])

# Fit multivariable logistic regression models to investigate odds ratio for missingness
# in poor hrqol at wave 4 between those with and without poor hrqol at earlier time points
glm(m_poor_hrqol_w4 ~ poor_hrqol_w1 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)
glm(m_poor_hrqol_w4 ~ poor_hrqol_w2 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)
glm(m_poor_hrqol_w4 ~ poor_hrqol_w3 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    family = "binomial", data = LSAC)

# Remove objects from global environment
rm("logit_2_w1", "logit_2_w2", "logit_2_w3", "logit_w1", "logit_w2", "logit_w3", 
   "LSAC", "miss_pattern", "ana_vars")

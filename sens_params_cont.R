# sens_params_cont
#
# This script will estimate relationships between HRQoL at age 10-11 years and 
# earlier measures of HRQoL. It is used to inform the choice of sensitivity 
# parameter values for the analysis of the continuous outcome.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
load(file = "../data/LSAC_AD.dta")

# Define set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# Add missingness indicators
LSAC$m_hrqol_w1 <- ifelse(is.na(LSAC$hrqol_w1), 1, 0)
LSAC$m_hrqol_w2 <- ifelse(is.na(LSAC$hrqol_w2), 1, 0)
LSAC$m_hrqol_w3 <- ifelse(is.na(LSAC$hrqol_w3), 1, 0)

# Examine missing data patterns 
miss_pattern <- mice::md.pattern(LSAC[, c("hrqol_w1", "hrqol_w2", "hrqol_w3", "hrqol_w4")], 
                                 rotate.names = TRUE)

# HRQoL summaries overall
LSAC[, c("hrqol_w1", "m_hrqol_w1", "hrqol_w2", "m_hrqol_w2", "hrqol_w3", "m_hrqol_w3")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         missing_text = "missing")

# HRQoL summaries by missingness in HRQoL at wave 4
LSAC[, c("hrqol_w1", "m_hrqol_w1", "hrqol_w2", "m_hrqol_w2", 
         "hrqol_w3", "m_hrqol_w3", "hrqol_w4", "m_hrqol_w4")] %>% 
  gtsummary::tbl_summary(type = all_dichotomous() ~ "categorical",
                         statistic = list(all_continuous() ~ "{mean} ({sd})",
                                          all_dichotomous() ~ "{n} ({p})"),
                         by = "m_hrqol_w4",
                         missing_text = "missing") 

# Fit univariable linear regression models to investigate mean difference in 
# hrqol at earlier waves between those with and without missing hrqol at wave 4
lm_w1 <- lm(hrqol_w1 ~ m_hrqol_w4, data = LSAC)
lm_w2 <- lm(hrqol_w2 ~ m_hrqol_w4, data = LSAC)
lm_w3 <- lm(hrqol_w3 ~ m_hrqol_w4, data = LSAC)
coef(lm_w1)[2]
coef(lm_w2)[2]
coef(lm_w3)[2]

# Fit multivariable linear regression models to investigate mean difference in 
# hrqol at earlier waves between those with and without missing hrqol at wave 4
lm(hrqol_w1 ~ m_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
   data = LSAC)
lm(hrqol_w2 ~ m_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
   data = LSAC)
lm(hrqol_w3 ~ m_hrqol_w4 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
   data = LSAC)

# Fit univariable logistic regression models to investigate log odds ratio for 
# missingness in hrqol at wave 4 for a unit increase in hrqol at earlier wave
glm(m_hrqol_w4 ~ hrqol_w1, data = LSAC, family = binomial(link = "logit"))
glm(m_hrqol_w4 ~ hrqol_w2, data = LSAC, family = binomial(link = "logit"))
glm(m_hrqol_w4 ~ hrqol_w3, data = LSAC, family = binomial(link = "logit"))

# Fit multivariable logistic regression models to investigate log odds ratio for 
# missingness in hrqol at wave 4 for a unit increase in hrqol at earlier wave
glm(m_hrqol_w4 ~ hrqol_w1 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    data = LSAC, family = binomial(link = "logit"))
glm(m_hrqol_w4 ~ hrqol_w2 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    data = LSAC, family = binomial(link = "logit"))
glm(m_hrqol_w4 ~ hrqol_w3 + overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, 
    data = LSAC, family = binomial(link = "logit"))

# Remove objects from global environment
rm("LSAC", "miss_pattern", "ana_vars", "lm_w1", "lm_w2", "lm_w3")
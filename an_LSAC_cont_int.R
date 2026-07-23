# an_LSAC_cont_int
#
# This script will conduct the extended sensitivity analysis for the 
# continuous outcome presented in Supplementary Material.
# This analysis includes an interaction between exposure and missingness.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
load(file = "../data/LSAC_AD.dta")
load(file = "results/cont_results.RData")

# Define set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# Select sensitivity parameters --------

# Add missingness indicators
LSAC$m_hrqol_w1 <- ifelse(is.na(LSAC$hrqol_w1), 1, 0)
LSAC$m_hrqol_w2 <- ifelse(is.na(LSAC$hrqol_w2), 1, 0)
LSAC$m_hrqol_w3 <- ifelse(is.na(LSAC$hrqol_w3), 1, 0)

# Fit univariable linear regression models to investigate mean difference in 
# hrqol at earlier waves between those with and without missing hrqol at wave 4
summary(lm(hrqol_w1 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4, data = LSAC))
summary(lm(hrqol_w2 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4, data = LSAC))
summary(lm(hrqol_w3 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4, data = LSAC))

# Fit multivariable linear regression models to investigate mean difference in 
# hrqol at earlier waves between those with and without missing hrqol at wave 4
summary(lm(hrqol_w1 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4 + sex + 
             indstat + noneng + age_w1 + sep_w1, data = LSAC))
summary(lm(hrqol_w2 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4 + sex + 
             indstat + noneng + age_w1 + sep_w1, data = LSAC))
summary(lm(hrqol_w3 ~ m_hrqol_w4 + overwt_w1 + overwt_w1:m_hrqol_w4 + sex + 
             indstat + noneng + age_w1 + sep_w1, data = LSAC))

# ------------------------------------------------------------------------------
# Delta-adjusted MI
# ------------------------------------------------------------------------------

cat("Running delta-adjusted MI")

# Sensitivity parameters
delta_df <- data.frame(delta1 = c("0.5", "-0.5", "-1.5"), 
                       delta2 = c("-2.5", "-1.5", "+0.5"))

# Set up meth and pred
meth <- c(rep("", 6), "mnar.norm")
names(meth) <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", 
                 "hrqol_w4")
pred <- rep(1, 7) %*% t(rep(1, 7)) - diag(7)
rownames(pred) <- colnames(pred) <- c("sex", "indstat", "noneng", "age_w1", 
                                      "sep_w1", "overwt_w1", "hrqol_w4")

# Estimate for each value of delta
for(i in 1:3){
  
  cat(".", "\n")
  
  imp <- mice(LSAC[, c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", 
                "hrqol_w4")],
       m = 50, maxit = 1, method = meth, pred = pred, 
       blots = list(hrqol_w4 = list(ums = paste(delta_df[i,1], delta_df[i,2], "*overwt_w1Overweight", sep = ""))), 
       seed = 70624, print = FALSE)
  
  # Mean HRQoL
  est <- summary(pool(with(imp, lm(hrqol_w4 ~ 1))), conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "mean",
                          meth = "Shifting",
                          param = paste("\u03b4\u2081", " = ", delta_df[i,1], 
                                        ", \u03b4\U2082 = ", delta_df[i, 2], sep = ""),
                          est = est$estimate,
                          ll = est$'2.5 %',
                          ul = est$'97.5 %'))
  
  # Regression-adjusted association of overwt with HRQoL
  est <- summary(pool(with(imp, lm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1))), conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "mean_diff",
                          meth = "Shifting",
                          param = paste("\u03b4\u2081", " = ", delta_df[i,1], 
                                        ", \u03b4\U2082 = ", delta_df[i, 2], sep = ""),
                          est =  est$estimate[2],
                          ll = est$'2.5 %'[2],
                          ul = est$'97.5 %'[2]))
}

cat("DONE", "\n")

# ------------------------------------------------------------------------------
# Extreme case analysis
# ------------------------------------------------------------------------------

# Extreme case: All missing HRQoL = min(HRQoL) overwt_w1 = "Overweight"; 100 otherwise
LSAC$hrqol_w4_ec <- ifelse(is.na(LSAC$hrqol_w4) & LSAC$overwt_w1 == "Overweight", 
                           min(LSAC$hrqol_w4, na.rm = TRUE),
                           ifelse(is.na(LSAC$hrqol_w4) & LSAC$overwt_w1 == "Not overweight",
                                  100, LSAC$hrqol_w4))

est <- lm(hrqol_w4_ec ~ 1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean",
                        meth = "Extreme case",
                        param = "miss=18.75 if overweight; else 100",
                        est = as.numeric(coef(est)),
                        ll = as.numeric(confint(est)[1]),
                        ul = as.numeric(confint(est)[2])))

est <- lm(hrqol_w4_ec ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean_diff",
                        meth = "Extreme case",
                        param = "miss=18.75 if overweight; else 100",
                        est = as.numeric(coef(est)[2]),
                        ll = as.numeric(confint(est)[2, 1]),
                        ul = as.numeric(confint(est)[2, 2])))

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

save("res", file = "results/cont_results_int.RData")

# ------------------------------------------------------------------------------
# Clean and graph results
# ------------------------------------------------------------------------------

# Filter by estimand
res_mean <- dplyr::filter(res, estimand == "mean")
res_mean_diff <- dplyr::filter(res, estimand == "mean_diff")

# Convert meth to an ordered factor and add exta rows for graphing
res_mean$meth <- factor(res_mean$meth, 
                        levels = c("Standard MI",
                                   "Shifting",
                                   "Weighting",
                                   "Extreme case"),
                        ordered = TRUE)
res_mean$meth <- fct_recode(res_mean$meth,
                            "Delta-adjusted MI" = "Shifting",
                            "Stacked MI" = "Weighting")
res_mean <- rbind(res_mean,
                  data.frame(estimand = rep(c("mean"), 3),
                             meth = c("**Primary analysis**",
                                      "**Sensitivity analysis**",
                                      "**Extensions**"),
                             param = rep("NA", 3),
                             est = rep(NA, 3),
                             ll = rep(NA, 3),
                             ul = rep(NA, 3)))

res_mean_diff$meth <- factor(res_mean_diff$meth, 
                             levels = c("Complete cases",
                                        "Shifting",
                                        "Weighting",
                                        "Extreme case"),
                             ordered = TRUE)
res_mean_diff$meth <- fct_recode(res_mean_diff$meth,
                                 "Delta-adjusted MI" = "Shifting",
                                 "Stacked MI" = "Weighting")

res_mean_diff <- rbind(res_mean_diff, 
                       data.frame(estimand = rep("mean_diff", 3),
                                  meth = c("**Primary analysis**",
                                           "**Sensitivity analysis**",
                                           "**Extensions**"),
                                  param = rep("NA", 3),
                                  est = rep(NA, 3),
                                  ll = rep(NA, 3),
                                  ul = rep(NA, 3)))

# Set up colours for graph
my_cols <- brewer.pal(n = 4, name = "Set1")

# Mean HRQoL ---------- 
p1 <- ggplot(res_mean) + 
  geom_point(aes(x = param, y = est, shape = meth, color = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = label_wrap(width = 10),
                   limits = c("\u03b4\u2081 = 0.5, \u03b4\U2082 = -2.5",
                              "\u03b4\u2081 = -0.5, \u03b4\U2082 = -1.5",
                              "\u03b4\u2081 = -1.5, \u03b4\U2082 = +0.5",
                              "**Extensions**",
                              "miss=18.75 if overweight; else 100",
                              "missing = 100",
                              "missing = 18.75",
                              "\u03C6 = -0.05",
                              "\u03C6 = -0.03",
                              "\u03C6 = -0.01",
                              "\u03C6 = -0.005",
                              "\u03b4 = -10",
                              "\u03b4 = -5",
                              "\u03b4 = -3",
                              "\u03b4 = -1",
                              "**Sensitivity analysis**",
                              "Standard MI",
                              "**Primary analysis**")) +
  scale_shape_manual(values = c("Standard MI" = 1,
                                "Delta-adjusted MI" = 17, 
                                "Stacked MI" = 15,
                                "Extreme case" = 16),
                     breaks = c("Standard MI",
                                "Delta-adjusted MI",
                                "Stacked MI",
                                "Extreme case")) +
  scale_colour_manual(values = c("Standard MI" = my_cols[1],
                                 "Stacked MI" = my_cols[2],
                                 "Delta-adjusted MI" = my_cols[3],
                                 "Extreme case" = my_cols[4]),
                      na.translate = FALSE) +
  labs(colour = "Approach", shape = "Approach") + 
  coord_flip() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 3, 1), "lines"),
        axis.text.y = element_markdown(),
        text = element_text(family = "serif")) +
  xlab(" ") + ylab("Mean HRQoL") +
  guides(shape = guide_legend(title = "Approach", 
                              override.aes = list(shape = c(1, 17, 15, 16))))
# Subscripts don't work in regular ggplot font. Change to "serif".
suppressWarnings(print(p1))
suppressWarnings(ggsave("results/mean_ext.jpg", plot = p1, width = 6.5, height = 4.7, unit = "in"))

# Association of being overweight with HRQoL ---------
p2 <- ggplot(res_mean_diff) + 
  geom_point(aes(x = param, y = est, shape = meth, colour = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10),
                   limits = c("\u03b4\u2081 = 0.5, \u03b4\U2082 = -2.5",
                              "\u03b4\u2081 = -0.5, \u03b4\U2082 = -1.5",
                              "\u03b4\u2081 = -1.5, \u03b4\U2082 = +0.5",
                              "**Extensions**",
                              "miss=18.75 if overweight; else 100",
                              "missing = 100",
                              "missing = 18.75",
                              "\u03C6 = -0.05",
                              "\u03C6 = -0.03",
                              "\u03C6 = -0.01",
                              "\u03C6 = -0.005",
                              "\u03b4 = -10",
                              "\u03b4 = -5",
                              "\u03b4 = -3",
                              "\u03b4 = -1",
                              "**Sensitivity analysis**",
                              "Complete cases",
                              "**Primary analysis**")) +
  scale_shape_manual(values = c("Complete cases" = 1,
                                "Delta-adjusted MI" = 17, 
                                "Stacked MI" = 15,
                                "Extreme case" = 16),
                     breaks = c("Complete cases",
                                "Delta-adjusted MI",
                                "Stacked MI",
                                "Extreme case")) +
  scale_colour_manual(values = c("Complete cases" = my_cols[1],
                                 "Delta-adjusted MI" = my_cols[2],
                                 "Stacked MI" = my_cols[3],
                                 "Extreme case" = my_cols[4]),
                      na.translate = FALSE) +
  labs(colour = "Approach", shape = "Approach") + 
  coord_flip() +
  theme_bw() +
  xlab(" ") + ylab("Difference in mean HRQoL") +
  guides(shape = guide_legend(title = "Approach", 
                              override.aes = list(shape = c(1, 17, 15, 16),
                                                  linetype = 1))) + 
  theme(axis.text.y = element_markdown(),
        plot.margin = unit(c(1, 1, 3, 1), "lines"),
        text = element_text(family = "serif"))
# Subscripts don't work in regular ggplot font. Change to "serif".
suppressWarnings(print(p2))
suppressWarnings(ggsave("results/mean_diff_ext.jpg", plot = p2, width = 6.5, height = 4.7, unit = "in"))

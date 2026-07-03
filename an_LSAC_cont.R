# an_LSAC_cont
#
# This script will conduct the sensitivity analysis for the continuous outcome.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
load(file = "../data/LSAC_AD.dta")

# Define set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# ------------------------------------------------------------------------------
# Primary analysis
# ------------------------------------------------------------------------------

# Set up meth and pred for mice
meth <- c(rep("", 6), "norm")
names(meth) <- c(ana_vars, "hrqol_w4")
pred <- rep(1, 7) %*% t(rep(1, 7)) - diag(7)
rownames(pred) <- colnames(pred) <- c(ana_vars, "hrqol_w4")
imp <- mice(LSAC[, c(ana_vars, "hrqol_w4")], maxit = 1, m = 50, seed = 70624,
            method = meth, predictorMatrix = pred, print = FALSE)

# Mean HRQoL using CCA: 77.6 (77.1, 78.0)
est <- lm(hrqol_w4 ~ 1, data = LSAC)
confint(est)
res <- data.frame(estimand = "mean", 
                  meth = "Complete cases",
                  param = "Complete cases",
                  est = as.numeric(coef(est)),                    
                  ll = confint(est)[1], 
                  ul = confint(est)[2])

# Mean HRQoL using MI: 77.2 (76.7, 77.6)
est <- pool(with(imp, lm(hrqol_w4 ~ 1)))
sum <- summary(est, conf.int = TRUE)
res <- data.frame(estimand = "mean",
                  meth = "Standard MI", 
                  param = "Standard MI",
                  est = sum$estimate,
                  ll = sum$'2.5 %', 
                  ul = sum$'97.5 %')

# Regression-adjusted association of overwt with HRQoL using CCA: -2.4 (-3.5, -1.4)
est <- lm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean_diff",
                        meth = "Complete cases", 
                        param = "Complete cases",
                        est = as.numeric(coef(est)[2]),
                        ll = confint(est)[2, 1],
                        ul = confint(est)[2, 2]))

# Regression-adjusted association of overwt with HRQoL using MI: -2.4 (-3.5, -1.4)
est <- pool(with(imp, lm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1)))
sum <- summary(est, conf.int = TRUE)
res <- rbind(res,
             data.frame(estimand = "mean_diff",
                        meth = "Standard MI", 
                        param = "Standard MI",
                        est = sum$estimate[2],
                        ll = sum$'2.5 %'[2],
                        ul = sum$'97.5 %'[2]))

# ------------------------------------------------------------------------------
# Delta-adjusted MI
# ------------------------------------------------------------------------------

cat("Running delta-adjusted MI")

# Sensitivity parameters
delta_vec <- c(-10, -5, -3, -1)

# Set up meth and pred
meth <- c(rep("", 6), "mnar.norm")
names(meth) <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", "hrqol_w4")
pred <- rep(1, 7) %*% t(rep(1, 7)) - diag(7)
rownames(pred) <- colnames(pred) <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", "hrqol_w4")

# Estimate for each value of delta
for(i in 1:length(delta_vec)){
  
  cat(".")
  
  mnar.blot <- list(hrqol_w4 = list(ums = paste(delta_vec[i])))
  imp <- mice(LSAC[, c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", "hrqol_w4")],
              m = 50, maxit = 1, method = meth, pred = pred, blots = mnar.blot, seed = 70624, print = FALSE)
  
  # Mean HRQoL
  est <- summary(pool(with(imp, lm(hrqol_w4 ~ 1))), conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "mean",
                          meth = "Shifting",
                          param = paste("\u03b4", " = ", delta_vec[i], sep = ""),
                          est = est$estimate,
                          ll = est$'2.5 %',
                          ul = est$'97.5 %'))
  
  # Regression-adjusted association of overwt with HRQoL
  est <- summary(pool(with(imp, lm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1))), conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "mean_diff",
                          meth = "Shifting",
                          param = paste("\u03b4", " = ", delta_vec[i], sep = ""),
                          est =  est$estimate[2],
                          ll = est$'2.5 %'[2],
                          ul = est$'97.5 %'[2]))
}

cat("DONE", "\n")

# ------------------------------------------------------------------------------
# Stacked MI
# ------------------------------------------------------------------------------

cat("Running stacked MI")

# Sensitivity parameters
phi_vec <- c(-0.005, -0.01, -0.03, -0.05)

# Set up meth and pred
meth <- c(rep("", 6), "norm", "")
names(meth) <- c(ana_vars, "hrqol_w4", "m_hrqol_w4")
pred <- rep(1, 8) %*% t(rep(1, 8)) - diag(8)
pred[8, ] <- pred[, 8] <- 0
rownames(pred) <- colnames(pred) <- c(ana_vars, "hrqol_w4", "m_hrqol_w4")

# Multiply impute and stack
imp <- mice(LSAC[, c(ana_vars, "hrqol_w4", "m_hrqol_w4")],
            m = 50, maxit = 1, method = meth, pred = pred, seed = 70624, print = FALSE)
stack <- complete(imp, action = "long", include = FALSE)

for(i in 1:length(phi_vec)){
  
  cat(".")
  
  # Calculate weights
  stack$wt <- ifelse(stack$m_hrqol_w4 == 1, exp(phi_vec[i] * stack$hrqol_w4), 1)
  stack <- as.data.frame(stack %>% group_by(.id) %>% dplyr::mutate(wt = wt / sum(wt)))
  
  # Check that the sum of weights over imputed data sets for each individual = 1
  # stack %>% group_by(.id) %>% dplyr::summarise(total_value = sum(wt))
  
  # Mean HRQoL
  fit <- glm(hrqol_w4 ~ 1, data = stack, family = gaussian(), weights = wt)
  coef(fit)
  jackcovar <- Jackknife_Variance(fit, stack, M = 50)
  VARIANCE_jack = diag(jackcovar)
  ci <- coef(fit) + c(-1, 1) * qnorm(0.975) * sqrt(VARIANCE_jack)
  res <- rbind(res,
               data.frame(estimand = "mean",
                          meth = "Weighting", 
                          param = paste("\u03C6", " = ", phi_vec[i], sep = ""),
                          est = as.numeric(coef(fit)),
                          ll = ci[1],
                          ul = ci[2]))
  
  # Regression-adjusted association of overwt with HRQoL
  fit <- glm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1,
             data = stack, family = gaussian(), weights = wt)
  coef(fit)[2]
  jackcovar <- StackImpute::Jackknife_Variance(fit, stack, M = 50)
  VARIANCE_jack = diag(jackcovar)
  ci <- coef(fit)[2] + c(-1, 1) * qnorm(0.975) * sqrt(VARIANCE_jack[2])
  res <- rbind(res,
               data.frame(estimand = "mean_diff",
                          meth = "Weighting", 
                          param = paste("\u03C6", " = ", phi_vec[i], sep = ""),
                          est = as.numeric(coef(fit)[2]),
                          ll = ci[1],
                          ul = ci[2]))
  
}

cat("DONE", "\n")

# ------------------------------------------------------------------------------
# Extreme case analysis
# ------------------------------------------------------------------------------

# Max and min values of HRQoL
summary(LSAC$hrqol_w4)

# Mean HRQoL ----

# Extreme case: All missing HRQoL = min(HRQoL): 
LSAC$hrqol_w4_ec1 <- ifelse(is.na(LSAC$hrqol_w4), 18.75, LSAC$hrqol_w4)
est <- lm(hrqol_w4_ec1 ~ 1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean",
                        meth = "Extreme case",
                        param = "missing = 18.75",
                        est = as.numeric(coef(est)),
                        ll = as.numeric(confint(est)[1]),
                        ul = as.numeric(confint(est)[2])))

# Extreme case: All missing HRQoL = 100: 
LSAC$hrqol_w4_ec2 <- ifelse(is.na(LSAC$hrqol_w4), 100, LSAC$hrqol_w4)
est <- lm(hrqol_w4_ec2 ~ 1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean",
                        meth = "Extreme case",
                        param = "missing = 100",
                        est = as.numeric(coef(est)),
                        ll = as.numeric(confint(est)[1]),
                        ul = as.numeric(confint(est)[2])))

# Association of being overweight with HRQoL ---- 

# Extreme case: All missing HRQoL = min(HRQoL)
est <- lm(hrqol_w4_ec1 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean_diff",
                        meth = "Extreme case",
                        param = "missing = 18.75",
                        est = as.numeric(coef(est)[2]),
                        ll = as.numeric(confint(est)[2, 1]),
                        ul = as.numeric(confint(est)[2, 2])))


# Extreme case: All missing HRQoL = 100
est <- lm(hrqol_w4_ec2 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1, data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "mean_diff",
                        meth = "Extreme case",
                        param = "missing = 100",
                        est = as.numeric(coef(est)[2]),
                        ll = as.numeric(confint(est)[2, 1]),
                        ul = as.numeric(confint(est)[2, 2])))

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

save("res", file = "results/cont_results.RData")

# ------------------------------------------------------------------------------
# Clean and graph results
# ------------------------------------------------------------------------------

load(file = "results/cont_results.RData")

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
res_mean <- rbind(res_mean,
                  data.frame(estimand = rep(c("mean"), 2),
                             meth = c("**Primary analysis**",
                                      "**Sensitivity analysis**"),
                             param = rep("NA", 2),
                             est = rep(NA, 2),
                             ll = rep(NA, 2),
                             ul = rep(NA, 2)))

res_mean_diff$meth <- factor(res_mean_diff$meth, 
                             levels = c("Complete cases",
                                        "Shifting",
                                        "Weighting",
                                        "Extreme case"),
                             ordered = TRUE)

res_mean_diff <- rbind(res_mean_diff, 
                       data.frame(estimand = rep("mean_diff", 2),
                                  meth = c("**Primary analysis**",
                                           "**Sensitivity analysis**"),
                                  param = rep("NA", 2),
                                  est = rep(NA, 2),
                                  ll = rep(NA, 2),
                                  ul = rep(NA, 2)))

# Set up colours for graph
my_cols <- brewer.pal(n = 4, name = "Set1")

# Mean HRQoL ---------- 
p1 <- ggplot(res_mean) + 
  geom_point(aes(x = param, y = est, shape = meth, color = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20),
                   limits = c("missing = 100",
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
                                "Shifting" = 17, 
                                "Weighting" = 15,
                                "Extreme case" = 16),
                     breaks = c("Standard MI",
                                "Shifting",
                                "Weighting",
                                "Extreme case")) +
  scale_colour_manual(values = c("Standard MI" = my_cols[1],
                                 "Shifting" = my_cols[2],
                                 "Weighting" = my_cols[3],
                                 "Extreme case" = my_cols[4]),
                      na.translate = FALSE) +
  labs(colour = "Approach", shape = "Approach") + 
  coord_flip() +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 3, 1), "lines"),
        axis.text.y = element_markdown()) +
  xlab(" ") + ylab("Mean HRQoL") +
  guides(shape = guide_legend(title = "Approach", 
                              override.aes = list(shape = c(1, 17, 15, 16))))
suppressWarnings(print(p1))
suppressWarnings(ggsave("mean.jpg", plot = p1, width = 6.5, height = 4.7, unit = "in"))

# Association of being overweight with HRQoL ---------
p2 <- ggplot(res_mean_diff) + 
  geom_point(aes(x = param, y = est, shape = meth, colour = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20),
                   limits = c("missing = 100",
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
                                "Shifting" = 17, 
                                "Weighting" = 15,
                                "Extreme case" = 16),
                     breaks = c("Complete cases",
                                "Shifting",
                                "Weighting",
                                "Extreme case")) +
  scale_colour_manual(values = c("Complete cases" = my_cols[1],
                                 "Shifting" = my_cols[2],
                                 "Weighting" = my_cols[3],
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
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
suppressWarnings(print(p2))
suppressWarnings(ggsave("results/mean_diff.jpg", plot = p2, width = 6.5, height = 4.7, unit = "in"))

# Remove objects from global environment 
rm("est", "fit", "imp", "LSAC", "jackcovar", "mnar.blot", "p1", "p2", "pred", 
   "res", "res_mean", "res_mean_diff", "stack", "sum", "ana_vars", "ci", 
   "delta_vec", "i", "meth", "my_cols", "phi_vec", "VARIANCE_jack")
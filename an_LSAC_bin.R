# an_LSAC_bin
#
# This script will conduct the sensitivity analysis for the binary outcome.
# 
# Written by R Mainzer
# ------------------------------------------------------------------------------

# Load data
load(file = "../data/LSAC_AD.dta")

# Define set of analysis variables excluding outcomes
ana_vars <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1")

# Recode poor_hrqol_w4 to have values 0 and 1
LSAC$poor_hrqol_w4 <- factor(1*(LSAC$poor_hrqol_w4 == "HRQoL < 50"))

# Number of multiple imputations
m <- 50

# Number of bootstrap replications
R <- 500

# Number of clusters for stacked MI (= number of individuals in data set)
n_clusters <- dim(LSAC)[1]

# ------------------------------------------------------------------------------
# Primary analysis
# ------------------------------------------------------------------------------

# Set up meth and pred for mice
meth <- c(rep("", 6), "logreg")
names(meth) <- c(ana_vars, "poor_hrqol_w4")
pred <- rep(1, 7) %*% t(rep(1, 7)) - diag(7)
rownames(pred) <- colnames(pred) <- c(ana_vars, "poor_hrqol_w4")
imp <- mice(LSAC[, c(ana_vars, "poor_hrqol_w4")], maxit = 1, m = m, 
            seed = 70624, method = meth, predictorMatrix = pred)

# Proportion with poor HRQoL using CCA: 5.4% (4.8%, 6.2%)
est <- glm(poor_hrqol_w4 ~ 1, family = binomial(link = "identity"), data = LSAC)
confint(est)
res <- data.frame(estimand = "prop",
                  meth = "Complete cases",
                  param = "Complete cases",
                  est = as.numeric(coef(est)),
                  ll = as.numeric(confint(est)[1]),
                  ul = as.numeric(confint(est)[2]))

# Proportion with poor HRQoL using MI: 6.0% (5.2%, 6.9%)
est <- pool(with(imp, glm(poor_hrqol_w4 ~ 1, 
                          family = binomial(link = "identity"))))
sum <- summary(est, conf.int = TRUE)
res <- rbind(res,
             data.frame(estimand = "prop",
                        meth = "Standard MI",
                        param = "Standard MI",
                        est = sum$estimate,
                        ll = sum$`2.5 %`,
                        ul = sum$`97.5 %`))

# Regression-adjusted association of overwt with poor HRQoL using CCA: 1.7% (-0.2%, 3.6%)
gcomp_cca <- gcomp_fun(data = LSAC, out = "poor_hrqol_w4", R = R)
res <- rbind(res,
             data.frame(estimand = "prop_diff",
                        meth = "Complete cases",
                        param = "Complete cases",
                        est = gcomp_cca$est,
                        ll = gcomp_cca$ll,
                        ul = gcomp_cca$ul))

# Regression-adjusted association of overwt with poor HRQoL using MI: 1.8% (-0.2%, 3.9%)
gcomp_mi <- gcomp_mi_fun(imp, out = "poor_hrqol_w4", m = m, R = R)
res <- rbind(res,
             data.frame(estimand = "prop_diff",
                        meth = "Standard MI",
                        param = "Standard MI",
                        est = gcomp_mi$est,
                        ll = gcomp_mi$ll,
                        ul = gcomp_mi$ul))

# ------------------------------------------------------------------------------
# Delta-adjusted MI
# ------------------------------------------------------------------------------

cat("Running delta-adjusted MI", "\n")

# Sensitivity parameters
delta_vec <- c(0.2, 0.5, 1, 2)

# Set up meth and pred
meth <- c(rep("", 6), "mnar.logreg")
names(meth) <- c("sex", "indstat", "noneng", "age_w1", "sep_w1", "overwt_w1", 
                 "poor_hrqol_w4")
pred <- rep(1, 7) %*% t(rep(1, 7)) - diag(7)
rownames(pred) <- colnames(pred) <- c("sex", "indstat", "noneng", "age_w1", 
                                      "sep_w1", "overwt_w1", "poor_hrqol_w4")

# Estimate for each value of delta
for(i in 1:length(delta_vec)){
  
  cat("delta = ", delta_vec[i], "\n")
  
  mnar.blot <- list(poor_hrqol_w4 = list(ums = paste(delta_vec[i])))
  imp <- mice(LSAC[, c("sex", "indstat", "noneng", "age_w1", "sep_w1", 
                       "overwt_w1", "poor_hrqol_w4")],
              m = m, maxit = 1, method = meth, pred = pred, 
              blots = mnar.blot, seed = 70624, print = FALSE)
  
  # Proportion with poor HRQoL
  est <- summary(pool(with(imp, glm(poor_hrqol_w4 ~ 1, 
                                    family = binomial(link = "identity")))), 
                 conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "prop",
                          meth = "Shifting",
                          param = paste("\u03b4", " = ", delta_vec[i], sep = ""),
                          est = est$estimate,
                          ll = est$'2.5 %',
                          ul = est$'97.5 %'))
  
  # Regression-adjusted association of overwt with poor HRQoL
  est <- gcomp_mi_fun(imp, out = "poor_hrqol_w4", m = m, R = R)
  #est <- summary(pool(with(imp, lm(hrqol_w4 ~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1))), conf.int = TRUE)
  res <- rbind(res,
               data.frame(estimand = "prop_diff",
                          meth = "Shifting",
                          param = paste("\u03b4", " = ", delta_vec[i], sep = ""),
                          est =  est[1],
                          ll = est[2],
                          ul = est[3]))
}

cat("DONE", "\n")

# ------------------------------------------------------------------------------
# Stacked MI
# ------------------------------------------------------------------------------

cat("Running stacked MI", "\n")

# Sensitivity parameters
phi_vec <- c(0.2, 0.5, 1, 2)

# Set up meth and pred
meth <- c(rep("", 6), "logreg", "")
names(meth) <- c(ana_vars, "poor_hrqol_w4", "m_hrqol_w4")
pred <- rep(1, 8) %*% t(rep(1, 8)) - diag(8)
pred[8, ] <- pred[, 8] <- 0
rownames(pred) <- colnames(pred) <- c(ana_vars, "poor_hrqol_w4", "m_hrqol_w4")

# Multiply impute and stack
imp <- mice(LSAC[, c(ana_vars, "poor_hrqol_w4", "m_hrqol_w4")],
            m = m, maxit = 1, method = meth, pred = pred, seed = 70624, print = FALSE)
stack <- complete(imp, action = "long", include = FALSE)

for(i in 1:length(phi_vec)){
  
  cat("phi = ", phi_vec[i])
  
  # Calculate weights
  stack$wt <- ifelse(stack$m_hrqol_w4 == 1, 
                     exp(phi_vec[i] * 1*(stack$poor_hrqol_w4 == 1)), 1)
  stack <- as.data.frame(stack %>% group_by(.id) %>% dplyr::mutate(wt = wt / sum(wt)))
  
  # Check that the sum of weights over imputed data sets for each individual = 1
  # stack %>% group_by(.id) %>% dplyr::summarise(total_value = sum(wt))
  
  # Proportion with poor HRQoL
  # Produces warning due to use of weights
  fit <- glm(poor_hrqol_w4 ~ 1, data = stack, family = binomial(link = "identity"), 
             weights = wt)
  coef(fit)
  jackcovar <- Jackknife_Variance(fit, stack, M = m)
  VARIANCE_jack = diag(jackcovar)
  ci <- coef(fit) + c(-1, 1) * qnorm(0.975) * sqrt(VARIANCE_jack)
  res <- rbind(res,
               data.frame(estimand = "prop",
                          meth = "Weighting", 
                          param = paste("\u03C6", " = ", phi_vec[i], sep = ""),
                          est = as.numeric(coef(fit)),
                          ll = ci[1],
                          ul = ci[2]))
  
  # Regression-adjusted association of overwt with HRQoL
  b_ests <- boot_est(n_clusters = n_clusters, data = stack, out = "poor_hrqol_w4", R = R)
  hist(b_ests$res) # Check distribution of bootstrap estimates
  dim(b_ests$b) # Check last bootstrap data set has correct dimension (243900 by 10)
  sd.boot <- sd(b_ests$res)
  mean.boot <- mean(b_ests$res)
  ci <- mean.boot + c(-1, 1) * 1.96 * sd.boot
  
  # Regression-adjusted association of overwt with HRQoL
  res <- rbind(res,
               data.frame(estimand = "prop_diff",
                          meth = "Weighting", 
                          param = paste("\u03C6", " = ", phi_vec[i], sep = ""),
                          est = mean.boot,
                          ll = ci[1],
                          ul = ci[2]))
  
  cat("\n")
  
}

cat("DONE", "\n")

# ------------------------------------------------------------------------------
# Extreme case analysis
# ------------------------------------------------------------------------------

# Extreme case: All missing HRQoL = poor HRQoL
LSAC$poor_hrqol_w4_ec1 <- factor(ifelse(is.na(LSAC$poor_hrqol_w4), 
                                       1, as.character(LSAC$poor_hrqol_w4)))

# Proportion with poor HRQoL
est <- glm(poor_hrqol_w4_ec1 ~ 1, family = binomial(link = "identity"), data = LSAC)
confint(est)
res <- rbind(res,
             data.frame(estimand = "prop",
                        meth = "Extreme case",
                        param = "missing = poor HRQoL",
                        est = as.numeric(coef(est)),
                        ll = as.numeric(confint(est)[1]),
                        ul = as.numeric(confint(est)[2])))

# Regression-adjusted association of overwt with HRQoL
est <- gcomp_fun(data = LSAC[c("poor_hrqol_w4_ec1", ana_vars)], 
                 out = "poor_hrqol_w4_ec1", R = R)
res <- rbind(res,
             data.frame(estimand = "prop_diff",
                        meth = "Extreme case",
                        param = "missing = poor HRQoL",
                        est = est$est,
                        ll = est$ll,
                        ul = est$ul))

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

save("res", file = "results/bin_results.RData")

# ------------------------------------------------------------------------------
# Clean and graph results
# ------------------------------------------------------------------------------

load(file = "results/bin_results.RData")

# Filter by estimand
res_prop <- dplyr::filter(res, estimand == "prop")
res_prop_diff <- dplyr::filter(res, estimand == "prop_diff")

# Convert meth to an ordered factor and add exta rows for graphing
res_prop$meth <- factor(res_prop$meth, 
                        levels = c("Standard MI",
                                   "Shifting",
                                   "Weighting",
                                   "Extreme case"),
                        ordered = TRUE)
res_prop <- rbind(res_prop,
                  data.frame(estimand = rep(c("prop"), 2),
                             meth = c("**Primary analysis**",
                                      "**Sensitivity analysis**"),
                             param = rep("NA", 2),
                             est = rep(NA, 2),
                             ll = rep(NA, 2),
                             ul = rep(NA, 2)))

res_prop_diff$meth <- factor(res_prop_diff$meth, 
                             levels = c("Complete cases",
                                        "Shifting",
                                        "Weighting",
                                        "Extreme case"),
                             ordered = TRUE)

res_prop_diff <- rbind(res_prop_diff, 
                       data.frame(estimand = rep("prop_diff", 2),
                                  meth = c("**Primary analysis**",
                                           "**Sensitivity analysis**"),
                                  param = rep("NA", 2),
                                  est = rep(NA, 2),
                                  ll = rep(NA, 2),
                                  ul = rep(NA, 2)))

# Set up colours for graph
my_cols <- brewer.pal(n = 4, name = "Set1")

# Proportion with poor HRQoL -----
p1 <- ggplot(res_prop) + 
  geom_point(aes(x = param, y = est, shape = meth, color = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20),
                   limits = c("missing = poor HRQoL",
                              "\u03C6 = 2",
                              "\u03C6 = 1",
                              "\u03C6 = 0.5",
                              "\u03C6 = 0.2",
                              "\u03b4 = 2",
                              "\u03b4 = 1",
                              "\u03b4 = 0.5",
                              "\u03b4 = 0.2",
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
  xlab("") + ylab("Proportion with poor HRQoL") +
  guides(shape = guide_legend(title = "Approach", 
                              override.aes = list(shape = c(1, 17, 15, 16))))
suppressWarnings(print(p1))
suppressWarnings(ggsave("results/prop.jpg", plot = p1, width = 6.5, height = 4.7, unit = "in"))

# Association of being overweight with HRQoL ---------
p2 <- ggplot(res_prop_diff) + 
  geom_point(aes(x = param, y = est, shape = meth, colour = meth)) +
  geom_linerange(aes(x = param, y = est, ymin = ll, ymax = ul, colour = meth),
                 key_glyph = "path") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20),
                   limits = c("missing = poor HRQoL",
                              "\u03C6 = 2",
                              "\u03C6 = 1",
                              "\u03C6 = 0.5",
                              "\u03C6 = 0.2",
                              "\u03b4 = 2",
                              "\u03b4 = 1",
                              "\u03b4 = 0.5",
                              "\u03b4 = 0.2",
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
  xlab("") + ylab("Difference in proportions of poor HRQoL") +
  guides(shape = guide_legend(title = "Approach", 
                              override.aes = list(shape = c(1, 17, 15, 16),
                                                  linetype = 1))) + 
  theme(axis.text.y = element_markdown(),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))
suppressWarnings(print(p2))
suppressWarnings(ggsave("results/prop_diff.jpg", plot = p2, width = 6.5, height = 4.7, unit = "in"))

# Remove objects from global environment 
rm("est", "fit", "imp", "LSAC", "jackcovar", "mnar.blot", "p1", "p2", "pred", 
   "res", "res_prop", "res_prop_diff", "stack", "sum", "ana_vars", "ci", 
   "delta_vec", "i", "meth", "my_cols", "phi_vec", "VARIANCE_jack")


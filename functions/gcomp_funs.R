# gcomp_funs
# Functions used to do g-computation

gcomp <- function(data, ind, out, glm.weights){
  # Estimate difference in risks using gcomp
  # data: data frame
  # ind: required 2nd argument for bootstrap function
  # out: outcome (string)
  
  data <- data[ind,]
  
  fit <- glm(as.formula(paste(out, "~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1")),
             family = binomial(link = "logit"), weights = glm.weights, data = data)
  
  mean_1 <- mean(predict(fit, 
                         newdata = replace(data, "overwt_w1", as.factor("Overweight")), 
                         type = "response"))
  mean_2 <- mean(predict(fit, 
                         newdata = replace(data, "overwt_w1", as.factor("Not overweight")),
                         type = "response"))
  mean_1 - mean_2
}

gcomp_fun <- function(data, out, R, glm.weights = NULL){
  # Use bootstrap with gcomp function to obtain SE and confidence interval
  
  bstrap <- boot(data = data, statistic = gcomp, stype = "i", out = out, glm.weights = glm.weights, R = R)
  
  bstrap_ci <- boot.ci(bstrap, conf = 0.95, type = "perc")
  
  data.frame(cbind(est = bstrap$t0, 
                   set = sqrt(var(bstrap$t)[1,1]),
                   ll = bstrap_ci$percent[4],
                   ul = bstrap_ci$percent[5]))
  
  #data.frame(cbind(est = bstrap$t0, 
  #                 se = sqrt(var(bstrap$t)[1,1]),
  #                 ll = with(bstrap, 2*t0 - mean(t) + qnorm(0.025) * sqrt(var(bstrap$t)[1,1])),
  #                 ul = with(bstrap, 2*t0 - mean(t) + qnorm(0.975) * sqrt(var(bstrap$t)[1,1]))))
  #boot.ci(bstrap, type = "norm")
  
}

gcomp_mi_fun <- function(imp, out, m, R){
  g.res <- lapply(c(1:m), 
                  function(i){
                    gcomp_fun(complete(imp, i), out = out, R = R)
                    }) %>%
    ldply(., data.frame)
  est <- mean(g.res$est)
  b <- var(g.res$est)
  v <- mean(g.res$se^2)
  t <- v + (1 + 1 / m) * b
  df <- (m - 1) * (1 + v / (b + b / m))^2
  data.frame(cbind(est = est,
                   ll = est - qt(0.975, df)*sqrt(t),
                   ul = est + qt(0.975, df)*sqrt(t)))
}
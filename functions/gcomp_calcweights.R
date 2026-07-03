# Function used to do g-comp that includes the calculation of the weights
# Note, weight calculation can be done outside this function for the cluster 
# bootstrap
gcomp_calcweights <- function(data, ind, out){
  # Estimate difference in risks using gcomp
  # data: data frame
  # ind: required 2nd argument for bootstrap function
  # out: outcome (string)
  
  data <- data[ind,]
  
  # Calculate weights
  data$wt <- ifelse(data$m_hrqol_w4 == 1, 
                     exp(phi_vec * 1*(data$poor_hrqol_w4 == 1)), 1)
  data <- as.data.frame(data %>% group_by(.id) %>% dplyr::mutate(wt = wt / sum(wt)))
  
  fit <- glm(as.formula(paste(out, "~ overwt_w1 + sex + indstat + noneng + age_w1 + sep_w1")),
             family = binomial(link = "logit"), weights = wt, data = data)
  
  mean_1 <- mean(predict(fit, 
                         newdata = replace(data, "overwt_w1", as.factor("Overweight")), 
                         type = "response"))
  mean_2 <- mean(predict(fit, 
                         newdata = replace(data, "overwt_w1", as.factor("Not overweight")),
                         type = "response"))
  mean_1 - mean_2
}
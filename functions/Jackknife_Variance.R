# RM modified version of L Beesley Jackknife_Variance function

Jackknife_Variance = function(fit, stack, M){
  if('glm' %in% class(fit)){
    if(substr(fit$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
      dispersion = 1
    }else{
      dispersion = StackImpute::glm.weighted.dispersion(fit)
    }
    covariance_weighted = summary(fit)$cov.unscaled*dispersion
  }else{
    covariance_weighted = vcov(fit)
  }
  results <- apply(cbind(c(1:M)), 1,FUN=StackImpute::func.jack, stack)
  #Nobs = length(stack[,1])
  #results_corrected = matrix(rep(as.vector(coef(fit)),M), ncol = M, byrow=F) - ((Nobs-M)/Nobs)*results
  if(is.matrix(results) == FALSE){
    theta_var = var(results)*(M-1)*((M-1)/M)
  } else {
    theta_var = var(t(results))*(M-1)*((M-1)/M)
  }
  Variance =covariance_weighted + (1+M)*theta_var
  return(Variance)
}
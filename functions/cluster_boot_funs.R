# cluster_boot_funs
# Functions used to do the cluster bootstrap

# Function to create one bootstrap sample
cr_boot_sample <- function(n_clusters, data){
  cluster_list <- split(data, data$.id)
  boot_ids <- sample(1:n_clusters, n_clusters, replace = TRUE)
  resampled_list <- cluster_list[boot_ids]
  names(resampled_list) <- paste0("boot_", 1:n_clusters, "_cluster_", boot_ids)
  boot_sample <- do.call(rbind, resampled_list)
  return(boot_sample)
}

# Function to create R bootstrap samples and use g-comp to estimate the 
# difference in risk for each bootstrap sample
boot_est <- function(n_clusters, data, out, R){
  
  res <- rep(0, R)
  
  for(i in 1:R){
    cat(".")
    if(i %% 25 == 0){cat(i, "\n")}
    b <- cr_boot_sample(n_clusters, data = stack)
    res[i] <- gcomp_calcweights(data = b, out = "poor_hrqol_w4")
  }
  
  list(res = res, b = b)
  
}
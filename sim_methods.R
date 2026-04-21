

# function to perform nlmr analyses using simulated data
sim_methods <- function(dat, rep_num) {


  # Method 0: linear MR
  #===================================================
  
  # collect summary results for linear mr
  #------------------------------------------------
  # regress x on g
  model_x <- lm(X ~ G, data = dat)
  summary_x <- summary(model_x)
  bx <- summary_x$coefficients["G", "Estimate"]
  bxse <- summary_x$coefficients["G", "Std. Error"]
  
  # regress y on g
  model_y <- lm(Y ~ G, data = dat)
  summary_y <- summary(model_y)
  by <- summary_y$coefficients["G", "Estimate"]
  byse <- summary_y$coefficients["G", "Std. Error"]
  
  # get summary stats for x
  xmean <- mean(dat$X)
  xmin <- min(dat$X)
  xmax <- max(dat$X)
  
  # organise the summary stats data
  m0 <- data.frame(
    bx = bx,
    bxse = bxse,
    by = by,
    byse = byse,
    xmean = xmean,
    xmin = xmin,
    xmax = xmax
  ) %>% 
    mutate(lace = by / bx,
           lace_se = byse / bx,
           rep = rep_num,
           method = "lmr",
           strata = 0)
  
  
  # Method 1: Doubly ranked method with X
  #===================================================
  m1 <- create_nlmr_summary2025(
    y = dat$Y,
    x = dat$X,
    g = dat$G,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>% 
    mutate(lace = by / bx,
           lace_se = byse / bx,
           rep = rep_num,
           method = "rank+X",
           strata = c(1:10))
  
  
  # Method 2: Doubly ranked method with X-GV
  #===================================================
  m2 <- create_nlmr_summary2025(
    y = dat$Y,
    x = dat$X,
    xs = dat$X1,
    g = dat$G,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>% 
    mutate(lace = by / bx,
           lace_se = byse / bx,
           rep = rep_num,
           method = "rank+(X-GV)",
           strata = c(1:10))
  
  # combine all method results for this replication
  results_all <- do.call(rbind, list(m0, m1, m2))
  
  # return all results
  return(results_all)

}





# Modified version of create_nlmr_summary:
# - Adds `xs` parameter to allow stratification on a variable other than exposure `x` (defaults to `x` if NULL).
# - Includes strata assignments in the output list as `strata`.
# - Uses `dplyr` for strata computation, enhancing readability and efficiency.
# - Ensures consistent matrix handling with `drop = FALSE` in subsetting.


# modified function
#===============================================================================

create_nlmr_summary2025 <- function(y,
                                    x,
                                    g,
                                    covar = NULL,
                                    family = "gaussian",
                                    controlsonly = FALSE,
                                    q,
                                    prestrat = 1,
                                    strata_method = "ranked",
                                    strata_bound = c(0.2, 0.1, 0.8, 0.9),
                                    extra_statistics = FALSE,
                                    report_GR = FALSE,
                                    report_het = FALSE,
                                    seed = 1234,
                                    xs = NULL) {
  
  # Preserve random seed state
  if (exists(".Random.seed")) {
    old <- .Random.seed
    on.exit({ .Random.seed <<- old })
  }
  if (!is.na(seed)) { set.seed(seed) }
  
  # Input validation checks
  stopifnot("report_GR only works with strata_method ranked" = !(report_GR == TRUE & strata_method != "ranked"))
  stopifnot("y must be a Surv object with family coxph" = !(family == "coxph" & !survival::is.Surv(y)))
  stopifnot("cannot use controlsonly option with family coxph" = !(family == "coxph" & controlsonly == TRUE))
  
  # Handle covariates if provided
  if (!is.null(covar) & !(is.matrix(covar) & is.numeric(covar))) {
    warning("covariates should be entered as numeric matrix: attempting to convert", immediate. = TRUE)
    covar2 <- model.matrix(~., data = as.data.frame(covar), na.action = na.pass)[, -1]
    print(head(covar2))
    user_input <- readline("Do you want to run using this matrix for covariates (y/n) ")
    if (user_input != 'y') stop('Exiting since you did not press y')
    covar <- as.matrix(covar2)
  }
  
  # Calculate the strata
  if (family == "binomial" | family == "gaussian" | family == "coxph") {
    if (strata_method == "ranked") {
      # Step 1: Compute strata1 from g for GR_stats
      z = rank(g, ties.method = "random")
      strata1 = floor((z - 1) / q / prestrat) + 1
      GR_stats <- getGRvalues(X = x, Zstratum = strata1)
      
      # Step 2: Create x0q based on xs or x
      if (is.null(xs)) {
        rank_var <- x
      } else {
        rank_var <- xs
      }
      id = seq_along(rank_var)
      temp <- data.frame(rank_var = rank_var, strata1 = strata1, id = id)
      temp <- dplyr::arrange(temp, rank_var)
      temp <- dplyr::group_by(temp, strata1)
      temp <- dplyr::mutate(temp, x0q = ceiling(rank(rank_var, ties.method = "random") / prestrat))
      temp <- dplyr::arrange(temp, id)
      x0q <- temp$x0q
    } else if (strata_method == "residual") {
      family2 = ifelse(family == "binomial", "binomial", "gaussian")
      ivf <- iv_free(y = y, x = x, g = g, covar = covar, q = q, family = family2, controlsonly = controlsonly)
      x0q <- ivf$x0q
    } else {
      stop("strata_method must be 'ranked' or 'residual'")
    }
  } else {
    stop("family must be gaussian or binomial or coxph")
  }
  
  quant <- q
  
  # Calculate associations for each stratum
  by <- rep(NA, quant)
  byse <- rep(NA, quant)
  bx <- rep(NA, quant)
  bxse <- rep(NA, quant)
  xmean <- rep(NA, quant)
  xmax <- rep(NA, quant)
  xmin <- rep(NA, quant)
  strata_stats <- vector("list", quant)
  
  for (j in 1:quant) {
    # Describe the quantiles of original data
    if (j == 1) {
      xmin[j] <- quantile(x[x0q == j], strata_bound[1])
    } else {
      xmin[j] <- quantile(x[x0q == j], strata_bound[2])
    }
    if (j == quant) {
      xmax[j] <- quantile(x[x0q == j], strata_bound[3])
    } else {
      xmax[j] <- quantile(x[x0q == j], strata_bound[4])
    }
    xmean[j] <- mean(x[x0q == j])
    
    # Model the beta coefficients based on family
    if (family == "binomial") {
      if (is.null(covar)) {
        model <- glm(y[x0q == j] ~ g[x0q == j], family = "binomial")
        if (controlsonly == TRUE) {
          model2 <- lm(x[x0q == j & y == 0] ~ g[x0q == j & y == 0])
        } else {
          model2 <- lm(x[x0q == j] ~ g[x0q == j])
        }
      } else {
        model <- glm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE], family = "binomial")
        if (controlsonly == TRUE) {
          model2 <- lm(x[x0q == j & y == 0] ~ g[x0q == j & y == 0] + covar[x0q == j & y == 0, , drop = FALSE])
        } else {
          model2 <- lm(x[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE])
        }
      }
      if (is.na(model$coef[2])) {
        stop("the regression coefficient of the outcome on the instrument in one of the quantiles is missing")
      }
      by[j] <- model$coef[2]
      byse[j] <- summary(model)$coef[2, 2]
    } else if (family == "coxph") {
      if (is.null(covar)) {
        model <- survival::coxph(y[x0q == j] ~ g[x0q == j])
        model2 <- lm(x[x0q == j] ~ g[x0q == j])
      } else {
        model <- survival::coxph(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE])
        model2 <- lm(x[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE])
      }
      if (is.na(model$coef[1])) {
        stop("the regression coefficient of the outcome on the instrument in one of the quantiles is missing")
      }
      by[j] <- model$coef[1]
      byse[j] <- summary(model)$coef[1, 3]
    } else {
      if (is.null(covar)) {
        model <- lm(y[x0q == j] ~ g[x0q == j])
        model2 <- lm(x[x0q == j] ~ g[x0q == j])
      } else {
        model <- lm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE])
        model2 <- lm(x[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = FALSE])
      }
      if (is.na(model$coef[2])) {
        stop("the regression coefficient of the outcome on the instrument in one of the quantiles is missing")
      }
      by[j] <- model$coef[2]
      byse[j] <- summary(model)$coef[2, 2]
    }
    
    if (is.na(model2$coef[2])) {
      stop("the regression coefficient of the exposure on the instrument in one of the quantiles is missing")
    }
    
    bx[j] <- model2$coef[2]
    bxse[j] <- summary(model2)$coef[2, 2]
    
    # Extra statistics if requested
    if (extra_statistics) {
      stats <- list(
        strata = j,
        xmin = quantile(x[x0q == j], 0),
        xmax = quantile(x[x0q == j], 1),
        xmean = mean(x[x0q == j]),
        ymin = quantile(y[x0q == j], 0),
        ymax = quantile(y[x0q == j], 1),
        x_fstat = (bx[j] / bxse[j])^2
      )
      strata_stats[[j]] <- stats
    }
    
    model <- NULL
    model2 <- NULL
  }
  
  # Prepare output data
  output <- data.frame(bx, by, bxse, byse, xmean, xmin, xmax)
  names(output) <- c("bx", "by", "bxse", "byse", "xmean", "xmin", "xmax")
  
  final_output_list = list(summary = output, strata = x0q)
  
  # Add extra statistics to output if requested
  if (extra_statistics) {
    stats <- as.data.frame(do.call(rbind, strata_stats))
    final_output_list[["strata_statistics"]] <- stats
  }
  
  # Add ranked method-specific outputs
  if (strata_method == "ranked") {
    final_output_list[["GR_max"]] <- GR_stats[1]
    
    # Test of IV-exposure assumption
    xcoef_sub <- bx
    xcoef_sub_se <- bxse
    p_het <- 1 - pchisq(metafor::rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE, df = (q - 1))
    p_het_trend <- metafor::rma.uni(xcoef_sub ~ xmean, vi = xcoef_sub_se^2, method = "DL")$pval[2]
    
    if (report_GR == TRUE) {
      final_output_list[["GR_results"]] <- GR_stats
    }
    if (report_het == TRUE) {
      p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
      final_output_list[["Heterogeneity_results"]] <- p_heterogeneity
    }
  }
  
  invisible(final_output_list)
}



# original function
#===============================================================================

# create_nlmr_summary <- function(y,
#                                 x,
#                                 g,
#                                 covar = NULL,
#                                 family = "gaussian",
#                                 controlsonly=FALSE,
#                                 q,
#                                 prestrat=1,
#                                 strata_method="ranked",
#                                 strata_bound=c(0.2,0.1,0.8,0.9),
#                                 extra_statistics =FALSE,
#                                 report_GR=FALSE,
#                                 report_het=FALSE,
#                                 seed=1234) {
#   
#   if( exists(".Random.seed") ) {
#     old <- .Random.seed
#     on.exit( { .Random.seed <<- old } )
#   }
#   if (!is.na(seed)) { set.seed(seed) }
#   
#   
#   # checks
#   stopifnot(
#     "report_GR only works with strata_method ranked" = !(report_GR==TRUE &
#                                                            strata_method!="ranked")
#   )
#   # coxph checks
#   stopifnot(
#     "y must be a Surv object with family coxph" = !(family=="coxph" &
#                                                       !is.Surv(y))
#   )
#   stopifnot(
#     "cannot use controlsonly option with family coxph" = !(family=="coxph" &
#                                                              controlsonly==TRUE)
#   )
#   
#   # covar issue
#   if (!is.null(covar) & !(is.matrix(covar) & is.numeric(covar))) {
#     warning("covariates should be entered as numeric matrix:
#             attempting to covert", immediate.=TRUE)
#     covar2<-model.matrix(~.,data=as.data.frame(covar), na.action=na.pass)[,-1]
#     print(head(covar2))
#     user_input <- readline("Do you want to run using this matrix for covariates (y/n) ")
#     if(user_input != 'y') stop('Exiting since you did not press y')
#     covar<-as.matrix(covar2)
#     
#   }
#   
#   
#   # calculate the iv-free association
#   if (family=="binomial" |family=="gaussian"| family=="coxph") {
#     if (strata_method=="residual"){
#       family2= ifelse(family=="binomial", "binomial", "gaussian")
#       ivf <- iv_free(
#         y = y, x = x, g = g,
#         covar = covar, q = q, family = family2, controlsonly=controlsonly
#       )
#       x0q <- ivf$x0q
#     } else if(strata_method=="ranked") {
#       # haodong ranked strata method
#       z = rank(g, ties.method = "random")
#       strata1 = floor((z-1)/q/prestrat)+1
#       # check GR statistic
#       GR_stats<-getGRvalues(X=x, Zstratum=strata1)
#       
#       id= seq(x)
#       temp<- data.frame(x=x,strata1=strata1,id=id, g=g)
#       temp<- arrange(.data=temp, x)
#       temp<-group_by(.data=temp, strata1)
#       temp<-mutate(.data=temp, x0q= ceiling(rank(x, ties.method = "random")/prestrat))
#       temp<-arrange(.data=temp, id)
#       
#       x0q <- temp$x0q
#     } else {
#       stop("strata ordering must be ranked or residual")
#     }
#   } else {
#     stop("family must be gaussian or binomial or coxph")
#   }
#   
#   
#   quant <- q
#   
#   # this calculates the association for each quanta
#   by <- rep(NA, quant)
#   byse <- rep(NA, quant)
#   bx <- rep(NA, quant)
#   bxse <- rep(NA, quant)
#   xmean <- rep(NA, quant)
#   xmax <- rep(NA, quant)
#   xmin <- rep(NA, quant)
#   true_xmax <- rep(NA, quant)
#   true_xmin <- rep(NA, quant)
#   strata_stats<-vector("list", quant)
#   
#   #use the ivfree quantiles
#   
#   for (j in 1:quant) {
#     # describe the quantiles of original data
#     if (j==1){
#       xmin[j] <- quantile(x[x0q == j], strata_bound[1])
#     }else{
#       xmin[j] <- quantile(x[x0q == j], strata_bound[2])
#     }
#     if (j==quant){
#       xmax[j] <- quantile(x[x0q == j], strata_bound[3])
#     }else{
#       xmax[j] <- quantile(x[x0q == j], strata_bound[4])
#     }
#     xmean[j] <- mean(x[x0q == j])
#     # model the beta coefficients
#     if (family == "binomial") {
#       if (is.null(covar)) {
#         model <- glm(y[x0q == j] ~ g[x0q == j], family = "binomial")
#         if (controlsonly==T){
#           model2 <- lm(x[x0q == j& y==0] ~ g[x0q == j & y==0])
#         } else {
#           model2 <- lm(x[x0q == j] ~ g[x0q == j])
#         }
#       }else{
#         model <- glm(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F],
#                      family = "binomial")
#         if (controlsonly==T){
#           model2 <- lm(x[x0q == j& y==0] ~ g[x0q == j & y==0]+ covar[x0q == j & y==0, , drop = F])
#         } else {
#           model2 <- lm(x[x0q == j] ~ g[x0q == j, drop = F] + covar[x0q == j, , drop = F])
#         }
#       }
#       if (is.na(model$coef[2])) {
#         stop("the regression coefficient of the outcome on the instrument
#            in one of the quantiles is missing")
#       }
#       by[j] <- model$coef[2]
#       byse[j] <- summary(model)$coef[2, 2]
#     }else if(family == "coxph"){
#       if (is.null(covar)) {
#         model <- coxph(y[x0q == j] ~ g[x0q == j])
#         model2 <- lm(x[x0q == j] ~ g[x0q == j])
#       }else{
#         model <- coxph(y[x0q == j] ~ g[x0q == j] + covar[x0q == j, , drop = F])
#         model2 <- lm(x[x0q == j] ~ g[x0q == j, drop = F] + covar[x0q == j, , drop = F])
#       }
#       if (is.na(model$coef[1])) {
#         stop("the regression coefficient of the outcome on the instrument
#            in one of the quantiles is missing")
#       }
#       by[j] <- model$coef[1]
#       byse[j] <- summary(model)$coef[1, 3]
#       
#     }else {
#       if (is.null(covar)) {
#         model <- lm(y[x0q == j] ~ g[x0q == j])
#         model2 <- lm(x[x0q == j] ~ g[x0q == j])
#       }else{
#         model <- lm(y[x0q == j] ~ g[x0q == j]+   covar[x0q == j, , drop = F])
#         model2 <- lm(x[x0q == j] ~ g[x0q == j]+   covar[x0q == j, , drop = F])
#       }
#       if (is.na(model$coef[2])) {
#         stop("the regression coefficient of the outcome on the instrument
#            in one of the quantiles is missing")
#       }
#       by[j] <- model$coef[2]
#       byse[j] <- summary(model)$coef[2, 2]
#     }
#     
#     if (is.na(model2$coef[2])) {
#       stop("the regression coefficient of the exposure on the instrument
#            in one of the quantiles is missing")
#     }
#     
#     bx[j] <- model2$coef[2]
#     bxse[j] <- summary(model2)$coef[2, 2]
#     
#     
#     if (extra_statistics) {
#       stats<- list( strata=j,
#                     xmin = quantile(x[x0q == j], 0),
#                     xmax = quantile(x[x0q == j], 1),
#                     xmean = mean(x[x0q == j]),
#                     ymin = quantile(y[x0q == j], 0),
#                     ymax = quantile(y[x0q == j], 1),
#                     x_fstat = (bx[j]/bxse[j])^2
#                     
#       )
#       strata_stats[[j]]<-append(strata_stats[[j]], stats)
#       
#     }
#     
#     
#     model <- NULL
#     model2 <- NULL
#   }
#   # output data
#   
#   output <- data.frame(bx, by, bxse, byse, xmean, xmin, xmax)
#   names(output) <- c("bx", "by", "bxse", "byse", "xmean", "xmin", "xmax")
#   
#   final_output_list= list(summary=output)
#   if (extra_statistics) {
#     stats<- as.data.frame(do.call(rbind, strata_stats))
#     final_output_list[["strata_statistics"]]<- stats
#   }
#   if (strata_method=="ranked"){
#     final_output_list[["GR_max"]]<- GR_stats[1]
#     
#     ##### Test of IV-exposure assumption #####
#     xcoef_sub <- bx
#     xcoef_sub_se <- bxse
#     p_het <- 1 - pchisq(rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE,
#                         df = (q - 1)
#     )
#     p_het_trend <- rma.uni(xcoef_sub ~ xmean,
#                            vi = xcoef_sub_se^2,
#                            method = "DL"
#     )$pval[2]
#     
#     if (report_GR==TRUE){
#       final_output_list[["GR_results"]]<-GR_stats
#     }
#     if (report_het==TRUE){
#       p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
#       final_output_list[["Heterogeneity_results"]]<- p_heterogeneity
#     }
#   }
#   
#   
#   # print(list(summary = head(output)))
#   invisible(final_output_list)
#   
# }



# function to check variables

custom_summary <- function(var) {
  if (is.numeric(var)) {
    mean_val <- mean(var, na.rm = TRUE)
    median_val <- median(var, na.rm = TRUE)
    sd_val <- sd(var, na.rm = TRUE)
    n_missing <- sum(is.na(var))
    return(list(type = "numeric", stats = c(mean = mean_val, median = median_val, sd = sd_val, n_missing = n_missing)))
  } else {
    counts <- table(var, useNA = "always")
    percentages <- prop.table(counts) * 100
    summary_df <- data.frame(category = names(counts), count = as.vector(counts), percentage = as.vector(percentages))
    summary_df$category <- ifelse(is.na(summary_df$category), "Missing", summary_df$category)
    return(list(type = "categorical", summary = summary_df))
  }
}



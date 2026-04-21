
# function to simulate the data for the falsification test
sim_nulloutcome <- function(g, 
                            exp,
                            n = length(g), 
                            bux = -1,
                            buy = 1) {
  
  # Ensure g and exp are provided and have the same length
  if (length(g) != length(exp)) {
    stop("g and exp must have the same length")
  }
  
  # Generate u and ey using faux::rnorm_multi with correlation 0 (independent)
  generated_data <- faux::rnorm_multi(
    n = n, 
    vars = 2, 
    mu = c(2, 2), 
    sd = c(1, 1), 
    r = 0, 
    varnames = c("u", "ey")
  )
  
  # Extract u and ey from the generated data
  u <- generated_data$u
  ey <- generated_data$ey
  
  # Compute x as exp + bux * u
  x <- exp + bux * u
  
  # Compute y as buy * u + epsilon_y
  y <- buy * u + ey
  
  # Create a data frame with the results
  dat <- data.frame(g = g, x = x, y = y)
  
  return(dat)
  
}

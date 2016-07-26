## Useful functions for HBV modelling and analysis

# R. T. Gray

# This script contains functions useful functions for plotting
# the HBV model results and for figure generation. 


midyearIndex <- function(n, timestep) {
  # This function extracts the mid year index. This used for 
  # calculating per population number annual outputs by dividing
  # annual value by mid year population size
  # 
  # Args:
  #   n: length of vector to extract midyear points
  #   timestep: modelling timestep used
  # Returns:
  #   indices associated with the mid year value
  #   
  #-----------------------------------------------------------------------
  
  seq(1, n, 1 / timestep) + round(1 / timestep / 2)
}

sumYear <- function(vector, timestep) {
  # This function is used to sum every 1/timestep values of a vector. 
  # Primarily used to calculate annual indicators
  # 
  # Args:
  #   vector: vector to sum
  #   timestep: modelling timestep used
  # Returns:
  #   Vector of sums
  # 
  #-----------------------------------------------------------------------
  
  return(unname(tapply(vector, (seq_along(vector)-1) %/% (1/timestep),
                       sum)))
}

yearDf <- function(df, sumVar, years, npops, timestep, 
                   midVar = NULL, divideVar = NULL) {
  # This function applys sumYear and midyear Index to a data frame of 
  # modelling results. Primarily used to produce a data frame with annual 
  # indicators rather than for each time step. 
  # 
  # Args:
  #   df:
  #   sumVar: String naming variable to calculate yearly sums
  #   years: Vector of years simulated
  #   npops: Number of populations in the dataframe
  #   timestep: Simulation timestep
  #   midVar: String naming variable to extract mid year estimates from.
  #   Default is NULL as it is optional.
  #   divideVar: String naming variable given by sumVar / midVar. Default
  #     is NULL as it is optional. Requires midVar to be non-Null
  # Returns:
  #   Dataframe with annual/yearly estimates
  # 
  # ----------------------------------------------------------------------
  
  yearFrame <- df[seq(1, nrow(df), 1 / timestep), ]
  yearFrame[, sumVar] <- sumYear(as.matrix(df[, sumVar]), timestep)
  yearFrame$year <- rep(head(years, -1), npops)
  
  if (!is.null(midVar)) {
    yearFrame[, midVar] <- df[midyearIndex(nrow(df), timestep), midVar]
  }
  
  if (!is.null(divideVar)) {
    yearFrame[, divideVar] <- yearFrame[, sumVar] / 
      yearFrame[, midVar]
  }
  return(yearFrame)
}

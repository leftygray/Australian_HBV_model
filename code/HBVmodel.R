## R function that simulates the HBV model equations

# This function simulates the HBV model updating the populations and key 
# indicators each time step from an initial population using a difference
# equation version of the model equations. 

# Authors: Richard T. Gray, Neil Bretana

HBVmodel <- function(pg, pm, initialPop, pts, transitions, 
                     interactions = NULL) {
  # Function to simulate the HBV model equations 
  # 
  # Args:
  #   pg: list containing project specifications
  #   pm: data frame with model parameters across columns and the value
  #     at each timestep down the rows
  #   initialPop: A named matrix giving the initial populaton sizes
  #     for each population and state. Each row represents an different 
  #     population and each column represents a different state
  #   transitions: A named matrix showing the rates people move from 
  #     one population to another. The rows show the population people move
  #     from and the columns show the population people move to.
  #   interactions: A named matrix showing which populations have 
  #     at-risk interactions
  #   pts: The simulation time points in years 
  # Returns: 
  #   A list containing the following outputs:
  #     allPops: A 3D array containing the population size in each
  #      population and state for each time point
  #     newInfections: A matrix containing the number of new infections
  #      each time step for each population
  #     newHBVdeaths:  A matrix containing the number of HBV deaths each 
  #       time step for each population
  #     newVaccinations: A matrix containing the number of new vaccinations
  #       each time step for each population
  #     newMigrants: A matrix containing the number of new migrants
  #       entering the population each time step for each population
  #     newTreatments: A matrix containing the number of new treatments
  #       each time step for each population
  #     newCured: A matrix containing the number of people who clear their
  #       infection each time step for each population
  #
  # TODO: 
  #   newDiagnoses (??)
  #   Make internal calculations population agnostic
  #   Fix up FOI equations/calculations
  #   Still some hard coded numbers which should be moved into parameters
  #   
  # -----------------------------------------------------------------------
  
  # Extract useful inputs -------------------------------------------------
  
  dt <- pg$timestep
  npops <- pg$npops
  nstates <- pg$nstates
  years <- pts
  npts <- length(pts)
  
  stateNames <- pg$state_names
  popNames <- pg$population_names
  dimNames <- list(popNames, stateNames)
  
  # Initialize output -----------------------------------------------------
  
  # May results in allPops array - population, state, time
  allPops <- array(0, c(npops, nstates, npts), dimnames = dimNames) 
  allPops[, , 1] <- initialPop # the first sheet is the initial population

  initResultsMatrix <- matrix(0, nrow = npops, ncol = npts)
  rownames(initResultsMatrix) <- popNames
  
  newInfections <- initResultsMatrix
  newHBVdeaths <- initResultsMatrix
  newVaccinations <- initResultsMatrix
  newMigrants <- initResultsMatrix
  newTreatments <- initResultsMatrix
  newCured <- initResultsMatrix
  
  # Force of infection calculations ---------------------------------------
  # TODO: add FOI equations and parameters to model
  f_o_i_0 <- 0
  f_o_i_1 <- 0 
  f_o_i_2 <- 0
  f_o_i_3 <- 0
  
  forceInfection <- matrix(0, ncol = 1, nrow = npops)
  forceInfection[1] <- f_o_i_0
  forceInfection[2] <- f_o_i_1
  forceInfection[3] <- f_o_i_2
  forceInfection[4] <- f_o_i_3
  
  # Progression of HBV only affects acute and chronics --------------------
  progress <- matrix(0, ncol = npts + 1, nrow = npops)
  progress[1, ] <- pm$ac_res_rate * pm$prog_chron_0 #age group 0 to 4
  progress[2, ] <- pm$ac_res_rate * pm$prog_chron_1  #age group 5 to 14
  progress[3, ] <- pm$ac_res_rate * pm$prog_chron_2  #age group 15 to 44
  progress[4, ] <- pm$ac_res_rate * pm$prog_chron_3  #age group 45  
  
  # Births ----------------------------------------------------------------
  # births differ based on time but only applied to susceptibles and only 
  # to age group 0 so other age groups have to be 0
  births <- matrix(0, ncol = npts + 1, nrow = npops)
  births[1, ] <- pm$births
  
  # Mortality -------------------------------------------------------------
  
  # Mortality is split into background mortality and HBV mortality
  # bgMortality changes according to time, applied to all age groups but 
  # disregards HBV status
  bgMortality <-  matrix(0, ncol = npts + 1, nrow = npops)
  bgMortality[1, ] <- pm$bgmort0to4
  bgMortality[2, ] <- pm$bgmort5to14  
  bgMortality[3, ] <- pm$bgmort15to44
  bgMortality[4, ] <- pm$bgmort45
  
  # hbvMortality differs according to age group and HBV status acute and 
  # chronic only
  hbvMortality <- array(0, c(npops, nstates, npts +1), dimnames = dimNames)
  hbvMortality[1, "a", ] <- pm$ac_mort_rate_0
  hbvMortality[2, "a", ] <- pm$ac_mort_rate_1
  hbvMortality[3, "a", ] <- pm$ac_mort_rate_2
  hbvMortality[4, "a", ] <- pm$ac_mort_rate_3
  hbvMortality[1, "ch", ] <- pm$chr_mort_rate_0
  hbvMortality[2, "ch", ] <- pm$chr_mort_rate_1
  hbvMortality[3, "ch", ] <- pm$chr_mort_rate_2
  hbvMortality[4, "ch", ] <- pm$chr_mort_rate_3
  
  # Vaccination -----------------------------------------------------------
  # vac only applies to susceptibles and immune, differs based on age 
  # groups depending on vac_eff, vacc_prop, vacc_prog, and vacc_avail
  # TODO: move hard coded inputs into parmaters
  vac_avail <- numeric(npts + 1)
  vac_prog_0 <- numeric(npts + 1)
  vac_prog_1 <- numeric(npts + 1)
  
  vac_avail[years >= 1985] <- pm$vacc_tog[1]
  vac_prog_0[years >= 2000] <- 1
  vac_prog_1[(years >= 1998) & (years < 2012)] <- 1 
  
  vac <- matrix(0, ncol = npts + 1, nrow = npops)
  vac[0, ] <- pm$vac_eff_0 * pm$vac_prop_0 * vac_prog_0 *
    0 * vac_avail #age group 0 to 4
  vac[1, ] <- pm$vac_eff_1 * pm$vac_prop_1 * vac_prog_1 * 
    1 * vac_avail #age group 5 to 14
  vac[2, ] <- pm$vac_eff_2 * pm$vac_prop_2 * 
    2 * vac_avail #age group 15 to 44
  vac[3, ] <- pm$vac_eff_3 * pm$vac_prop_3 *
    3 * vac_avail #age group 45 
  
  # Migration -------------------------------------------------------------
  # migration proportions differ based on age group, HBV status, and time;
  # acute 0, mig_series is labeled migseries* in pm
  # mig_pred time? Where to put this line?
  mig_pred <- numeric(npts + 1)
  mig_pred[years < 2011] <- 1 / pm$mig_series[1]
  mig_pred[years >= 2011] <- 1

  migration <- array(0, c(npops, nstates, npts + 1), dimnames = dimNames) 
  migration["age0to4", "s", ] <- pm$mig0to4sus * pm$mig_series[1] *
    mig_pred
  migration["age5to14", "s", ] <- pm$mig5to14sus  * pm$mig_series[1] *
    mig_pred
  migration["age15to44", "s", ] <- pm$mig15to44sus * pm$mig_series[1] *
    mig_pred
  migration["age45", "s", ] <- pm$mig45sus * pm$mig_series[1] * mig_pred
  migration["age0to4", "ch", ] <- pm$mig0to4chronic * pm$mig_series[1] *
    mig_pred
  migration["age5to14", "ch", ] <- pm$mig5to14chronic * pm$mig_series[1] *
    mig_pred
  migration["age15to44", "ch", ] <- pm$mig15to44chronic * 
    pm$mig_series[1] * mig_pred
  migration["age45", "ch", ] <- pm$mig45chronic * pm$mig_series[1] *
    mig_pred
  migration["age0to4", "cl", ] <- pm$mig0to4cleared * pm$mig_series[1] *
    mig_pred
  migration["age5to14", "cl", ] <- pm$mig5to14cleared * pm$mig_series[1] *
    mig_pred
  migration["age15to44", "cl", ] <- pm$mig15to44cleared * 
    pm$mig_series[1] * mig_pred
  migration["age45", "cl", ] <- pm$mig45cleared * pm$mig_series[1] *
    mig_pred
  migration["age0to4", "i", ] <- 0 #* pm$mig_series[1] * mig_pred
  migration["age5to14", "i", ] <- 0 #* pm$mig_series[1] * mig_pred
  migration["age15to44", "i", ] <- 0 #* pm$mig_series[1] * mig_pred
  migration["age45", "i", ] <- 0 #* pm$mig_series[1] * mig_pred
    
  # Clearance -------------------------------------------------------------
  # clearance only applied to chronics and cleared, differ by age group
  clear <- matrix(0, nrow = npops, ncol = npts + 1)
  clear[1, ] <- pm$clr_rate_0 #age group 0 to 4
  clear[2, ] <- pm$clr_rate_1 #age group 5 to 14
  clear[3, ] <- pm$clr_rate_2 #age group 15 to 44
  clear[4, ] <- pm$clr_rate_3 #age group 45  
  
  # Recovery --------------------------------------------------------------
  recover <- matrix(0, nrow = npops, ncol = npts + 1)
  recover[1, ] <- pm$ac_res_rate * (1 - pm$prog_chron_0) 
  recover[2, ] <- pm$ac_res_rate * (1 - pm$prog_chron_1) 
  recover[3, ] <- pm$ac_res_rate * (1 - pm$prog_chron_2)  
  recover[4, ] <- pm$ac_res_rate * (1 - pm$prog_chron_3)  
  
  # Loop over state equations ---------------------------------------------
  
  for (time in 2:npts) {
    
    oldPop <- allPops[, , time - 1]
    newPop <- allPops[, , time]
    
    # Equations 

    newPop[, "s"] <- oldPop[, "s"] + 
      as.numeric(transitions %*% oldPop[, "s"]) +
      forceInfection * oldPop[, "s"] +
      births[, time] + migration[, "s", time] -
      bgMortality[, time] * oldPop[, "s"] -
      hbvMortality[, "s", time] * oldPop[, "s"] -
      vac[, time] * oldPop[, "s"]
    
    # For acutes there appears to be no migration. Probably makes sense
    # because acute stage is short so when people enter Australia they are 
    # either susceptible of chronically infected. May need review later. 
    newPop[, "a"] <-  oldPop[, "a"] + 
      as.numeric(transitions %*% oldPop[, "a"]) +
      forceInfection * oldPop[, "a"] +
      migration[, "a", time] -
      progress[, time] * oldPop[, "a"] -
      bgMortality[, time] * oldPop[, "a"] -
      hbvMortality[, "a", time] * oldPop[, "a"] -
      recover[, time] * oldPop[, "a"]         

    newPop[, "ch"] <- oldPop[, "ch"] + 
      as.numeric(transitions %*% oldPop[, "ch"]) +
      progress[, time] * oldPop[, "a"] +
      migration[, "ch", time] -
      bgMortality[, time] * oldPop[, "ch"] -
      hbvMortality[, "ch", time] * oldPop[, "ch"] -   
      clear[, time] * oldPop[, "ch"]
   
    newPop[, "cl"] <- oldPop[, "cl"] + 
      as.numeric(transitions %*% oldPop[, "cl"]) +
      migration[, "cl", time] -
      bgMortality[, time] * oldPop[, "cl"] -
      hbvMortality[, "cl", time] * oldPop[, "cl"] + 
      recover[, time] * oldPop[, "a"] + 
      clear[, time] * oldPop[, "ch"]     

    newPop[, "i"] <- oldPop[, "i"] + 
      as.numeric(transitions %*% oldPop[, "i"]) +
      migration[, "i", time] -
      bgMortality[, time] * oldPop[, "i"] -
      hbvMortality[, "i", time] * oldPop[, "i"] + 
      vac[, time] * oldPop[, "s"]   
    
    # Sort out results 
    allPops[, , time] <- newPop
    newInfections[, time] <- forceInfection * oldPop[, "s"]
    newHBVdeaths[, time] <- hbvMortality[, "a", time] * oldPop[, "a"] + 
      hbvMortality[, "ch", time] * oldPop[, "ch"]
    newVaccinations[, time] <- vac[, time] * oldPop[, "s"]
    newMigrants[, time] <- migration[, "s", time] + 
      migration[, "a", time] + 
      migration[, "ch", time] + migration[, "cl", time] + 
      migration[, "i", time]
    newTreatments[, time] <- clear[, time] * oldPop[, "ch"]
    newCured[, time] <- recover[, time] * oldPop[, "a"]  
  }
  
  # Sort out results ------------------------------------------------------
  
  results <- list(allPops = allPops, 
                  newInfections = newInfections, 
                  newHBVdeaths = newHBVdeaths, 
                  newVaccinations = newVaccinations, 
                  newMigrants = newMigrants, 
                  newTreatments = newTreatments, 
                  newCured = newCured)
  
  return(results)
  
}

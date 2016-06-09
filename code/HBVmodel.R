## R function that simulates the HBV model equations

# This function simulates the HBV model updating the populations and key 
# indicators each time step from an initial population using a difference
# equation version of the model equations. 

# Authors: Richard T. Gray, Neil Bratana

HBVmodel <- function(pg, pm, initialPop, transitions = NULL, 
                     interactions = NULL, pts) {
  # Function to simulate the HBV model equations 
  # 
  # Args:
  #   pg: list containing project specifications
  #   pm: data frame with model parameters across columns and the value
  #     at each timestep down the rows
  #   initialPop: A 2D matrix (named?) giving the initial populaton sizes
  #     for each population and state. Each row represents an different 
  #     population and each column represents a different state
  #   transitions: A 2D array (matrix?) showing the rates people move from 
  #     one population to another. The rows show the population people move
  #     from and the columns show the population people move to.
  #   interactions: A 2D array (matrix?) showing which populations have 
  #     at-risk interactions
  #   pts: Time points the model is simulated for (maybe unecessary)
  # Returns: 
  #   A list containing the following outputs:
  #     allPops: A 3D array showing the population size in each population
  #       and state for each time point
  #     newInfections:
  #     newHBVdeaths:
  #     newVaccinations: ??
  #     newMigrants: ??
  #     newTreatments:
  #     newCured:
  #
  # Todo: 
  #   newDiagnoses (??)
  #   
  # -----------------------------------------------------------------------
  
  # Extract useful inputs -------------------------------------------------
  
  dt <- pg$dt
  # npts <- length(pts)
  npops <- pg$npops
  nstates <- pg$nstates
  years <- pg$years # 1951 to 2050
  npts <- length(pts) # length(years)
  
  stateNames <- c("s", "a", "ch", "cl", "i")
  popNames <- c("age0to4", "age5to14", "age15to44", "age45")
  dimNames <- list(popNames, stateNames)
  
  # Initialize output -----------------------------------------------------
  
  # May results in allPops array - population, state, time
  allPops <- array(0, c(npops, nstates, npts), dimnames = dimNames) 
  allpops[, , 1] <- initialPop # the first sheet is the initial population

  initResultsMatrix <- matrix(0, nrow = npops, ncol = npts, dimnames = dimNames)
  
  newInfections <- initResultsMatrix
  newHBVdeaths <- initResultsMatrix
  newVaccinations <- initResultsMatrix
  newMigrants <- initResultsMatrix
  newTreatments <- initResultsMatrix
  newCured <- initResultsMatrix
  
  # Force of infection calculations ---------------------------------------
  forceInfection <- matrix(0, ncol = 1, nrow = npops)
  forceInfection[1] <- pm$f_o_i_0
  forceInfection[2] <- pm$f_o_i_1
  forceInfection[3] <- pm$f_o_i_2
  forceInfection[4] <- pm$f_o_i_3
  
  #progression of HBV only affects acute and chronics
  #prog <- matrix(0, ncol = 4, nrow = 4)
  progress <- matrix(0, ncol = 1, nrow = npops)
  progress[0] <- pm$ac_res_rate * pm$prog_chron_0 #age group 0 to 4
  progress[1] <- pm$ac_res_rate * pm$prog_chron_1  #age group 5 to 14
  progress[2] <- pm$ac_res_rate * pm$prog_chron_2  #age group 15 to 44
  progress[3] <- pm$ac_res_rate * pm$prog_chron_3  #age group 45  
  
  # births differ based on time but only applied to susceptibles and only 
  # to age group 0 so other age groups have to be 0
  births <- matrix(0, ncol = npops, nrow = npts)
  births[, 1] <- pm$births
  
  #Mortality is split into background mortality and HCV mortality
  #bgMortality changes according to time, applied to all age groups but 
  #disregards HBV status
  bgMortality <-  matrix(0, ncol = npops, nrow = npts)
  bgMortality[, 1] <- pm$bgmort0to4
  bgMortality[, 2] <- pm$bgmort5to14  
  bgMortality[, 3] <- pm$bgmort15to44
  bgMortality[, 4] <- pm$bgmort45
  
  #hbvMortality differs according to age group and HBV status acute and chronic only, disregards time
  hbvMortality <- array(0, c(npops, nstates), dimnames = dimNames)
  hbvMortality[1,"a"] <- pm$ac_mort_rate0
  hbvMortality[2,"a"] <- pm$ac_mort_rate1
  hbvMortality[3,"a"] <- pm$ac_mort_rate2
  hbvMortality[4,"a"] <- pm$ac_mort_rate3
  hbvMortality[1,"ch"] <- pm$chr_mort_rate0
  hbvMortality[2,"ch"] <- pm$chr_mort_rate1
  hbvMortality[3,"ch"] <- pm$chr_mort_rate2
  hbvMortality[4,"ch"] <- pm$chr_mort_rate3
  
  #vacc only applies to susceptibles and immune, differs based on age 
  #groups depending on vac_eff, vacc_prop, vacc_prog, and vacc_avail
  #vacc <- matrix(0, ncol = npops, nrow = 1)
  #vacc <- array(0, c(npops, nstates), dimnames = dimNames) 
  vacc <- matrix(0, ncol = 1, nrow = npops)
  vacc[0] <- pm$vacc_eff_0*pm$vacc_prop_0*pm$vacc_prog*0*pm$vacc_avail #age group 0 to 4
  vacc[1] <- pm$vacc_eff_1*pm$vacc_prop_1*pm$vacc_prog*1*pm$vacc_avail #age group 5 to 14
  vacc[2] <- pm$vacc_eff_2*pm$vacc_prop_2*pm$vacc_prog*2*pm$vacc_avail #age group 15 to 44
  vacc[3] <- pm$vacc_eff_3*pm$vacc_prop_3*pm$vacc_prog*3*pm$vacc_avail #age group 45 
  
  #migration proportions differ based on age group, HBV status, and time; acute 0
  #mig_series is labeled migseries* in pm
  #mig_pred time? Where to put this line?
  mig_pred <- matrix(1, 1, npts)
  mig_pred(years < 2011) <- 1/pm$mig_series

  migration <- array(0, c(npops, nstates, npts), dimnames = dimNames) 
  migration["age0to4", "s"] <- pm$mig0to4sus * pm$mig_series * pm$mig_pred
  migration["age5to14", "s"] <- pm$mig5to14sus  * pm$mig_series * pm$mig_pred
  migration["age15to44", "s"] <- pm$mig15to44sus * pm$mig_series * pm$mig_pred
  migration["age45", "s"] <- pm$mig45sus * pm$mig_series * pm$mig_pred
  migration["age0to4", "ch"] <- pm$mig0to4chronic * pm$mig_series * pm$mig_pred
  migration["age5to14", "ch"] <- pm$mig5to14chronic * pm$mig_series * pm$mig_pred
  migration["age15to44", "ch"] <- pm$mig15to44chronic * pm$mig_series * pm$mig_pred
  migration["age45", "ch"] <- pm$mig45chronic * pm$mig_series * pm$mig_pred
  migration["age0to4", "cl"] <- pm$mig0to4cleared * pm$mig_series * pm$mig_pred
  migration["age5to14", "cl"] <- pm$mig5to14cleared * pm$mig_series * pm$mig_pred
  migration["age15to44", "cl"] <- pm$mig15to44cleared * pm$mig_series * pm$mig_pred
  migration["age45", "cl"] <- pm$mig45cleared * pm$mig_series * pm$mig_pred
  migration["age0to4", "i"] <- 0 * pm$mig_series * pm$mig_pred
  migration["age5to14", "i"] <- 0 * pm$mig_series * pm$mig_pred
  migration["age15to44", "i"] <- 0 * pm$mig_series * pm$mig_pred
  migration["age45", "i"] <- 0 * pm$mig_series * pm$mig_pred
    
  #clearance only applied to chronics and cleared, differ by age group
  #clear <- matrix(pm$clr_rate_1, pm$clr_rate_2, pm$clr_rate_3, 
  #pm$clr_rate_4,nrow = npops, dimnames = dimNames)
  clear <- matrix(0, ncol = 1, nrow = npops)
  clear[0] <- pm$clr_rate_0 #age group 0 to 4
  clear[1] <- pm$clr_rate_1 #age group 5 to 14
  clear[2] <- pm$clr_rate_2 #age group 15 to 44
  clear[3] <- pm$clr_rate_3 #age group 45  
  
  #recover
  recover <- matrix(acc, ncol = 1, nrow = npops)
  recover[0] <- pm$ac_res_rate*(1-pm$prog_chron_0) #age group 0 to 4
  recover[1] <- pm$ac_res_rate*(1-pm$prog_chron_1) #age group 5 to 14
  recover[2] <- pm$ac_res_rate*(1-pm$prog_chron_2) #age group 15 to 44
  recover[3] <- pm$ac_res_rate*(1-pm$prog_chron_3) #age group 45  
  
  # Transitions - setup for moving between populations due to aging or
  # changing characteristics. WARNING: hard coded for the time being
  transitions <- matrix(0, npops, npops)
  transitions[1,2] <- 1/5
  transitions[2,3] <- 1/10
  transitions[3,4] <- 1/30
  
  # Loop over state equations ---------------------------------------------
  
  for (time in 2:npts) {
    
    oldPop <- allpops[, , time-1]
    newPop <- allpops[, , time]
    
    # Equations 
    #added + before forceInfection 
    newPop[, "s"] <- oldPop[, "s"] + 
                    transistions %*% oldPop[, "s"] +
                    forceInfection * oldPop[, "s"] +
                    births[time, ] + migration[, "s", time] -
                    bgMortality[time, ] -
                    hbvMortality[, "s"] -
                    vacc * oldPop[, "s"]
    # For acutes there appears to be no migration. Probably makes sense
    # because acute satge is short so when people enter Australia they are 
    # either susceptible of chronically infected. May need review later. 
    newPop[, "a"] <-  oldPop[, "a"] + 
      transistions * oldPop[, "a"] +
                    forceInfection * oldPop[, "s"] +
                    migration[, "a", time] +
                    prog * oldPop[, "a"] -
                    bgMortality[time, ] -
                    hbvMortality[, "a"] -
                    recover * oldPop[, "a"]         

    newPop[, "ch"] <- oldPop[, "ch"] + 
      transistions * oldPop[, "ch"] +
                    prog * oldPop[, "ch"] +
                    migration[, "ch", time] -
                    bgMortality[time, ] -
                    hbvMortality[, "ch"]  -   
                    clear * oldPop[, "ch"]
   
    newPop[, "cl"] <- oldPop[, "cl"] + 
      transistions * oldPop[, "cl"] +
                    migration[, "cl", time] -
                    bgMortality[time, ] -
                    hbvMortality[, "cl"] + 
                    recover * oldPop[, "a"] + 
                    clear * oldPop[, "ch"]     

    newPop[, "i"] <- oldPop[, "i"] + 
      transistions * oldPop[, "i"] +
                    migration[, "i", time] -
                    bgMortality[time, ] -
                    hbvMortality[, "i"] + 
                    vacc * oldPop[, "s"]   
    
    # Sort out results 
    allPops[, , time] <- newPop
    newInfections[, time] <- forceInfection * oldPop[, "s"]
    newHBVdeaths[, time] <- hbvMortality[, "a"] + hbvMortality[, "ch"] 
    newVaccinations[, time] <- vacc * oldPop[, "s"]
    newMigrants[, time] <- migration[, "s", time] + migration[, "a", time] + 
      migration[, "ch", time] + migration[, "cl", time] + 
      migration[, "i", time]
    newTreatments[, time] <- clear * oldPop[, "ch"]
    newCured[, time] <- recover * oldPop[, "a"]  
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

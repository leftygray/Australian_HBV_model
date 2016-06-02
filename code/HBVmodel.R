## R function that simulates the HBV model equations

# This function simulates the HBV model updating the populations and key 
# indicators each time step from an initial population using a difference
# equation version of the model equations. 

# Authors: Richard T. Gray, Neil Bratana

HbvModel <- function(pg, pm, initialPop, transitions = NULL, 
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
  npts <- length(pts)
  npops <- pg$npops
  nstates <- pg$nstates
  
  stateNames <- c("s", "a", "ch", "cl", "i")
  popNames <- c("age0to4", "age5to14", "age15to44", "age45")
  dimNames <- list(popNames, stateNames)
  
  # Initialize output -----------------------------------------------------
  
  # May results in allPops array - population, state, time
  allPops <- array(0, c(npops, nstates, npts), dimnames = dimNames) 
  allpops[, , 1] <- initialPop # the first sheet is the initial population

  initResultsMatrix <- matrix(0, nrow = npops, ncol = nstates, 
                              dimnames = dimNames)
  newInfections <- initResultsMatrix
  # newHBVdeaths <- initResultsMatrix
  # newTreatments <- initResultsMatrix
  # newCured <- initResultsMatrix
  
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
  progress[1] <- pm$clr_rate_1 * pm$prog_chron_1  #age group 5 to 14
  progress[2] <- pm$clr_rate_2 * pm$prog_chron_2  #age group 15 to 44
  progress[3] <- pm$clr_rate_3 * pm$prog_chron_3  #age group 45  
  
  # births differ based on time but only applied to susceptibles and only 
  # to age group 0 so other age groups have to be 0
  births <- matrix(0, ncol = npops, nrow = npts)
  births[, 1] <- pm$births
  
  #Mortality is split into background mortality and HCV mortality
  #bgMortality changes according to time, applied to all age groups but 
  #disregards HBV status
  bgMortality <-  matrix(0, ncol = npops, nrow = npts)
  
  #hbvMortality differs according to age group and HBV status, disregards time
  hbvMortality <- array(0, c(npops, nstates), dimnames = dimNames)
  
  #vacc only applies to susceptibles and immune, differs based on age 
  #groups depending on vac_eff, vacc_prop, vacc_prog, and vacc_avail
  #vacc <- matrix(0, ncol = npops, nrow = 1)
  #vacc <- array(0, c(npops, nstates), dimnames = dimNames) 
  vacc <- matrix(0, ncol = 1, nrow = 4)
  vacc[0] <- vacc_eff_0*vacc_prop_0*vacc_prog*0*vacc_avail #age group 0 to 4
  vacc[1] <- vacc_eff_1*vacc_prop_1*vacc_prog*1*vacc_avail #age group 5 to 14
  vacc[2] <- vacc_eff_2*vacc_prop_2*vacc_prog*2*vacc_avail #age group 15 to 44
  vacc[3] <- vacc_eff_3*vacc_prop_3*vacc_prog*3*vacc_avail #age group 45 
  
  #migration proportions differ based on age group, HBV status, and time
  migration <- array(0, c(npops, nstates, npts), dimnames = dimNames) 
  
  #clearance only applied to chronics and cleared, differ by age group
  #clear <- matrix(pm$clr_rate_1, pm$clr_rate_2, pm$clr_rate_3, 
  #pm$clr_rate_4,nrow = npops, dimnames = dimNames)
  clear <- matrix(0, ncol = 1, nrow = 4)
  clear[0] <- pm$clr_rate_0 #age group 0 to 4
  clear[1] <- pm$clr_rate_1 #age group 5 to 14
  clear[2] <- pm$clr_rate_2 #age group 15 to 44
  clear[3] <- pm$clr_rate_3 #age group 45  
  
  #recover
  recover <- matrix(acc, ncol = 1, nrow = 4)
  recover[0] <- ac_res_rate*(1-prog_chron_0) #age group 0 to 4
  recover[1] <- ac_res_rate*(1-prog_chron_1) #age group 5 to 14
  recover[2] <- ac_res_rate*(1-prog_chron_2) #age group 15 to 44
  recover[3] <- ac_res_rate*(1-prog_chron_3) #age group 45  
  
  
  # Loop over state equations ---------------------------------------------
  
  for (time in 2:npts) {
    
    oldPop <- allpops[, , time-1]
    newPop <- allpops[, , time]
    
    # Equations 
    #added + before forceInfection 
    newPop[, "s"] <- oldPop[, "s"] + 
                    transistion * oldPop[, "s"] +
                    forceInfection * oldPop[, "s"] +
                    births[time, ] + migration[, "s", time] -
                    bgMortality[time, ] -
                    hbvMortality[, "s"] -
                    vacc * oldPop[, "s"]
    # For acutes there appears to be no migration. Probably makes sense
    # because acute satge is short so when people enter Australia they are 
    # either susceptible of chronically infected. May need review later. 
    newPop[, "a"] <-  oldPop[, "a"] + 
                    transistion * oldPop[, "a"] +
                    forceInfection * oldPop[, "s"] +
                    migration[, "a", time] +
                    prog * oldPop[, "a"] -
                    bgMortality[time,] -
                    hbvMortality[, "a"] -
                    recover * oldPop[, "a"]         

    newPop[, "ch"] <- oldPop[, "ch"] + 
                    transistion * oldPop[, "ch"] +
                    prog * oldPop[, "ch"] +
                    migration[, "ch", time] -
                    bgMortality[time,] -
                    hbvMortality[, "ch"]  -   
                    clear * oldPop[, "ch"]
   
    newPop[, "cl"] <- oldPop[, "cl"] + 
                    transistion * oldPop[, "cl"] +
                    migration[, "cl", time] -
                    bgMortality[time,] -
                    hbvMortality[, "cl"] + 
                    recover * oldPop[, "a"] + 
                    clear * oldPop[, "ch"]     

    newPop[, "i"] <- oldPop[, "i"] + 
                    transistion * oldPop[, "i"] +
                    migration[, "i", time] -
                    bgMortality[time,] -
                    hbvMortality[, "i"] + 
                    vacc * oldPop[, "s"]   
    
    # Sort out results 
    allPops[, , time] <- newPop
    newInfections[, time] <- forceInfection * oldPop[, "s"]
    # chronics 
  }
  
  # Sort out results ------------------------------------------------------
  
  results <- list(allPops = allPops, 
                  newInfections = newInfections)
#  totalPopulation
#  susceptible
#  acute
  # chronic
  # cleared
  # immune
  
  return(results)
  
}

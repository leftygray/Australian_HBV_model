## Function used fo  calibrate the HBV model for a inputed project

# This function is used to load updated model paramters and rerun the HBV 
# model for the best estimate parameters. It's primary purpose is to 
# assist with the calibration of the model to available data. It is a
# combination of the 1-SetupModel.Rmd and 2-RunHBVmodel.Rmd scripts. 
#
# 1-SetupModel.Rmd needs to be run before using this function to ensure all # the project specifications are entered and after the calibrations is 
# complete to save the parmeters into the project.
# 
# Calibration process.....
#
# Authors: Richard T. Gray, Neil Bretana

CalibrateHIVmodel <- function(project, resource = FALSE) {
  # Function for calibrating a HBV model project 
  # 
  # Args:
  #   project: list containing project specifications
  #   resource: Set to true to source main libraries and functions. 
  #     Should only need to do this once.
  # Returns: 
  #   Generates a plot for visual inspection.
  #
  # -----------------------------------------------------------------------
  
  graphics.off() # Close current plot
  
  basePath <- getwd()
  project_directory <- file.path(basePath, "projects")
 
  
  if (resource) {
    Rcode <- file.path(basePath, "code")
    
    # Load useful libraries
    source(file.path(Rcode, "LoadLibrary.R"))
    source(file.path(Rcode, "DataLibraries.R"))
    source(file.path(Rcode, "PlotOptions.R"))
    
    # Source model files
    source(file.path(Rcode, "HBVmodel.R"))
  } 
  
  # Load current project specs and inputs ---------------------------------
  projectFolder <- file.path(project_directory, project)
  load(file.path(projectFolder, paste0(project_name, ".rda")))
  
  # Read in input files (which have been updated by the user) -------------
  initialPops <- read_csv(file.path(projectFolder, 
                                    "initial_populations.csv"))
  constants <- read_csv(file.path(projectFolder,
                                  "parameters_constants.csv"))
  timeVarying <- read_csv(file.path(projectFolder,
                                    "parameters_time_varying.csv"))
  timeFactors <- read_csv(file.path(projectFolder,
                                    "time_varying_ranges.csv"))
  transitions <- read_csv(file.path(projectFolder,
                                    "population_transitions.csv"))
  interactions <- read_csv(file.path(projectFolder,
                                     "population_interactions.csv"))
  
  # Convert transistions and interactions into matrices
  transitions <- as.matrix(transitions[, 2:(pg$npops + 1)])
  rownames(transitions) <- colnames(transitions)
  
  interactions <- as.matrix(interactions[, 2:(pg$npops + 1)])
  rownames(interactions) <- colnames(interactions)
  
  # Update best estimate parameter set ------------------------------------
  bestEstimateYears <- timeVarying
  
  # Reshape constants so they can be appended to timevarying
  constantsSpread <- constants %>%
    select(parameter, value)
  parameters <- constantsSpread$parameter
  constantsSpread <- as.data.frame(t(constantsSpread[,-1]))
  colnames(constantsSpread) <- parameters
  rownames(constantsSpread) <- NULL

  # Append to time varyingCC
  constantsDf <- constantsSpread[rep(1,nrow(timeVarying)),]
  bestEstimateYears <- cbind(timeVarying, constantsDf)
  
  # Expand to all points by linear extrapolation bewteen years
  best_estimates <- as.data.frame(matrix(0, ncol = ncol(bestEstimateYears),
                                         nrow = pg$npts))
  colnames(best_estimates) <- colnames(bestEstimateYears)
  
  for (var in colnames(bestEstimateYears)) {
    tempValues <- bestEstimateYears[, var]
    for (year in 1:(pg$nyears-1)) {
      indices <- ((year-1)* 1/pg$timestep + 1): (year * 1/pg$timestep + 1)
      yearValues <- seq(tempValues[year], tempValues[year+1], 
                        length = (1 + 1/pg$timestep))
      best_estimates[indices, var] <- yearValues
    }  
  }
  
  # Set up best estimate initial population
  init_pop <- filter(initialPops, parameter == "init_pop")$value
  
  initial_pop <- initialPops %>%
    filter(parameter != "init_pop")
  
  initialPopMat <- as.data.frame(matrix(initial_pop$value, 
                                        ncol = pg$nstates + 1,  
                                        nrow = pg$npops))
  colnames(initialPopMat) <- c(pg$state_names, "pop_prop")
  rownames(initialPopMat) <- pg$population_names
  
  popProp <- init_pop * initialPopMat$pop_prop
  
  best_initial_pop <- apply(initialPopMat[, 1:pg$nstates], 2, 
                            function(x) x * popProp)
  
  # Run model on updated best estimates -----------------------------------
  results <- HBVmodel(pg, best_estimates, best_initial_pop, pg$pts,
                          transitions)
  
  # Create calibration plot ----------------------------------------------
 
  # Plot in a separate figure
  windows(width = 35, height=30, xpos = 200) # Dimesnions big enough to 
                                             # maximize
  # Total population size 
  
  popSizes <- as.data.frame.table(bestResults$allPops)
  colnames(popSizes) <- c("age_group", "state", "time_step", "popsize") 
  popSizes$time_step <- as.numeric(popSizes$time_step)
  
  popSizes <- popSizes %>% # add pts as well
    mutate(year = pg$start_year + (time_step-1) * pg$timestep) 
  
  totalPop <- popSizes %>%
    group_by(year) %>%
    summarise(totalpop = sum(popsize)) 
  
  totalPopPlot <- ggplot(data = totalPop, aes(x = year, y = totalpop)) +
    geom_line(colour = "blue") +
    xlab("Year") + ylab("Population size") +
    scale_x_continuous(breaks = seq(pg$start_year, pg$end_year +1, 
                                    by = 20)) + plotOpts +
    ggtitle("Total population size")
  
  # Population size by age group
  popSizesAge <- popSizes %>%
    group_by(year, age_group) %>%
    summarise(totalpop = sum(popsize)) %>%
    ungroup() %>%
    arrange(age_group)

  plotAgePops <- popSizesAge
  plotAgePops$age_group <- factor(plotAgePops$age_group)
  
  # rename levels for plotting
  levels(plotAgePops$age_group) <- c("Age < 4", "Age 5 to 14", 
                                     "Age 15 to 44", "Age > 45") 
  
  agePopPlot <- ggplot(data = plotAgePops, 
                       aes(x = year, y = totalpop, group = age_group)) +
    geom_line(aes(colour = age_group)) + 
    scale_colour_brewer(name = "Age Group", palette = "Set1",
                        labels = c("Age < 4", "Age 5 to 14", 
                                   "Age 15 to 44",
                                   "Age > 45")) +
    scale_x_continuous(breaks = seq(pg$start_year, pg$end_year +1, 
                                    by = 50)) +
    xlab("Year") + ylab("Population size") +
    plotOpts + theme(legend.position = "right") + 
    ggtitle("Age group population size")
  
  # New Infections
  newInfections <- as.data.frame(t(bestResults$newInfections)) %>%
    mutate(year = pg$pts) %>%
    select(year, everything()) %>%
    gather("age_group", "infections", 2:(pg$npops+1)) %>%
    mutate(incidence = infections / popSizesAge$totalpop)
  
  newInfectionsPlot <- ggplot(data = newInfections, 
      aes(x = year, y = infections, group= age_group)) +
    geom_line(aes(color = age_group)) + 
    xlab("Year") + ylab("New infections") +
    scale_x_continuous(breaks = seq(pg$start_year, pg$end_year +1, 
                                    by = 20)) +
    scale_colour_brewer(name = "Age Group", palette = "Set1",
                        labels = c("Age < 4", "Age 5 to 14", 
                                   "Age 15 to 44",
                                   "Age > 45")) + 
    plotOpts + theme(legend.position = "right") +
    ggtitle("New infections by age")
    
  incidencePlot <- ggplot(data = newInfections, 
      aes(x = year, y = incidence * 1e5, group = age_group)) + 
    geom_line(aes(color = age_group)) + 
    xlab("Year") + ylab("Incidence per 100,000") +
    scale_x_continuous(breaks = seq(pg$start_year, pg$end_year +1, 
                                    by = 20)) +
    scale_colour_brewer(name = "Age Group", palette = "Set1",
                        labels = c("Age < 4", "Age 5 to 14", 
                                   "Age 15 to 44",
                                   "Age > 45")) +
    plotOpts + theme(legend.position = "right") +
    ggtitle("Incidence by age")
  
    # Put all plots in a grid
    print(plot_grid(totalPopPlot, agePopPlot, newInfectionsPlot, 
                    incidencePlot, ncol = 3, nrow = 3))
  
  
}

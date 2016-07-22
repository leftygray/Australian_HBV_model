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
  
  basePath <- getwd()
  project_directory <- file.path(basePath, "projects")
  
  if (resource) {
    Rcode <- file.path(basePath, "code")
    
    # Load useful libraries
    source(file.path(Rcode, "LoadLibrary.R"))
    source(file.path(Rcode, "DataLibraries.R"))
    
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
  
  # Following gives a warning but I think it is okay (maybe need
  # to convert constantsSpread into a repeated dataframe)
  bestEstimateYears <- cbind(timeVarying, constantsSpread[1, ])
  
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
  
  # Create simple plot ----------------------------------------------------
  
}

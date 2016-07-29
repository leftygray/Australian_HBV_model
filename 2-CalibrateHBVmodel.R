## Function used fo  calibrate the HBV model for a inputed project

# This function is used to load updated model paramters and rerun the HBV 
# model for the best estimate parameters. It's primary purpose is to 
# assist with the calibration of the model to available data. It is a
# combination of the 1-SetupModel.Rmd and 2-RunHBVmodel.Rmd scripts. 
#
# 1-SetupModel.Rmd needs to be run before using this function to ensure all
# the project specifications are entered and after the calibrations is 
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
  # Initilization --------------------------------------------------------
  
  graphics.off() # Close current plot
  
  # Need to set wd to source file location
  basePath <- getwd()
  project_directory <- file.path(basePath, "projects")
 
  
  if (resource) {
    Rcode <- file.path(basePath, "code")
    
    # Load useful libraries
    source(file.path(Rcode, "LoadLibrary.R"))
    source(file.path(Rcode, "DataLibraries.R"))
    source(file.path(Rcode, "PlotOptions.R"))
    source(file.path(Rcode, "PlotFunctions.R"))
    LoadLibrary(cowplot) # Note masks ggsave
    
    # Source model files
    source(file.path(Rcode, "HBVmodel.R"))
  } 
  
  #Load current project specs and inputs ---------------------------------
  projectFolder <- file.path(project_directory, project)
  load(file.path(projectFolder, paste0(project, ".rda")))
  
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
  waifw_matrix <- read_csv(file.path(projectFolder,
                                     "population_waifw_matrix.csv"))
  
  # Convert transistions and waifw_matrix into matrices
  transitions <- as.matrix(transitions[, 2:(pg$npops + 1)])
  rownames(transitions) <- colnames(transitions)
  
  waifw_matrix <- as.matrix(waifw_matrix[, 2:(pg$npops + 1)])
  rownames(waifw_matrix) <- colnames(waifw_matrix)
  
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
                          transitions, waifw_matrix)
  
  # Create calibration plot ----------------------------------------------
 
  # Plot in a separate figure
  windows(width = 35, height=30, xpos = 200) # Dimensions big enough to 
                                             # maximize
  # Population sizes
  totalPop <- popResults(pg, bestResults,
                         populations = "all", states = "all") 
  
  totalPopPlot <- indicatorPlot(totalPop, ylabel = "Population size",
                                range = FALSE, 
                                xlimits = c(pg$start_year, 
                                            pg$end_year, 20)) +
    ggtitle("Overall population")
  
  
  # Plot of population size in every group 
  popSizesAge <- popResults(pg, bestResults,
                            states = "all") %>%
    FactorPop(popLabels)
  
  agePopPlot <- indicatorPlot(popSizesAge, ylabel = "Population size",
                              xlimits = c(pg$start_year, 
                                          pg$end_year, 20),
                              range = FALSE, 
                              groupPlot = "population") +
    ggtitle("Population by age")
  
  # Prevalence 
  totalNumInfected <- popResults(pg, bestResults,
                                 states = c("a", "ch"), 
                                 populations = "all")
  
  ageNumInfected <- popResults(pg, bestResults, 
                               states = c("a", "ch")) %>%
    FactorPop(popLabels)
  
  totalInfectedPlot <- indicatorPlot(totalNumInfected, 
                                     ylabel = "Total number infected",
                                     range = FALSE, 
                                     xlimits = c(pg$start_year, 
                                                 pg$end_year, 20)) +
    ggtitle("Overall number infected")
  
  popInfectedPlot <- indicatorPlot(ageNumInfected, 
                                   ylabel = "Number infected",
                                   range = FALSE, 
                                   xlimits = c(pg$start_year, 
                                               pg$end_year, 20),
                                   groupPlot = "population") +
    ggtitle("Infected by population")
  
  # Total prevalance
  tempInfected <- popResults(pg, bestResults,
                             states = c("a", "ch"), 
                             range = FALSE, populations = "all") %>% 
    select(-year)
  
  tempTotal <- popResults(pg, bestResults, 
                          populations = "all", states = "all",
                          range = FALSE) %>% select(-year)
  
  tempPrev <- data.frame(year = totalPop$year, 
    best = (100 * tempInfected / tempTotal)) %>%
    tbl_df()
  
  totalPrevPlot <- indicatorPlot(tempPrev, 
                                 ylabel = "Overall prevalence (%)",
                                 range = FALSE, 
                                 xlimits = c(pg$start_year, 
                                             pg$end_year, 20)) +
    ggtitle("Overall prevalence")
  
  # Pop prevalence
  tempInfected <- popResults(pg, bestResults,
                             states = c("a", "ch"),
                             range = FALSE) %>% select(-year, -population)
  
  tempTotal <- popResults(pg, bestResults,
                          states = "all",
                          range = FALSE) %>% select(-year, -population)
  
  tempPrev <- data.frame(year = popSizesAge$year, 
                         population = popSizesAge$population,
                         best = (100 * tempInfected / tempTotal)) %>%
    tbl_df()
  
  popPrevPlot <- indicatorPlot(tempPrev, 
                               ylabel = "Population prevalence (%)",
                               xlimits = c(pg$start_year, 
                                           pg$end_year, 20),
                               range = FALSE, 
                               groupPlot = "population") +
    ggtitle("Prevalence by population")
  
  # New infections 
  
  # Annual number of new infections
  totalInfections <- indicatorResults(pg, bestResults,
                                      "newInfections",
                                      populations = "all", 
                                      annual = "sum")
  
  popInfections <- indicatorResults(pg, bestResults, 
                                    "newInfections",
                                    annual = "sum") %>%
    FactorPop(popLabels)
  
  totalNewInfectionsPlot <- indicatorPlot(totalInfections, 
                                     ylabel = "New infections",
                                     range = FALSE, 
                                     xlimits = c(pg$start_year, 
                                                 pg$end_year, 20)) +
    ggtitle("Overall new infections")
  
  popNewInfectionsPlot <- indicatorPlot(popInfections, 
                                        ylabel = "New infections",
                                        xlimits = c(pg$start_year, 
                                                    pg$end_year, 20),
                                        range = FALSE, 
                                        groupPlot = "population") +
    ggtitle("Overall new infections")
  
  # Cumulative incidence plots
  tempTotalInfections <- indicatorResults(pg, bestResults,
                                          "newInfections",
                                          populations = "all", 
                                          annual = "sum") 
  
  totalCumInfections <- tempTotalInfections %>%
    mutate(best = cumsum(best)) 
  
  totalCumInfectionsPlot <- indicatorPlot(totalCumInfections, 
                                      ylabel = "Number of infections",
                                      range = FALSE, 
                                      xlimits = c(pg$start_year, 
                                                  pg$end_year, 20)) +
    ggtitle("Overall cumulative infections")
  
  tempPopInfections <- indicatorResults(pg, bestResults,
                                        "newInfections", 
                                        annual = "sum") %>%
    FactorPop(popLabels)
  
  popCumInfections <- tempPopInfections %>%
    group_by(population) %>%
    mutate(best = cumsum(best)) 
  
  popCumInfectionsPlot <- indicatorPlot(popCumInfections, 
                                        ylabel = "Number of infections",
                                        range = FALSE, 
                                        xlimits = c(pg$start_year, 
                                                    pg$end_year, 20),
                                        groupPlot = "population") +
    ggtitle("Cumulative infections by  population")
  
  # Incidence 
  
  tempTotalInfections <- select(tempTotalInfections, -year)
  tempPopInfections <- select(tempPopInfections, -year)
  
  tempTotalPop <- popResults(pg, bestResults,
                             populations = "all", states = "all",
                             range = FALSE) %>% select(-year)
  tempTotalPop <- tempTotalPop[MidyearIndex(nrow(tempTotalPop),
                                            pg$timestep), ]
  
  
  tempPopPop <- popResults(pg, bestResults, states = "all",
                           range = FALSE) %>% 
    select(-year) %>%
    arrange(population) %>%
    FactorPop(popLabels)
  tempPopPop <- tempPopPop[MidyearIndex(nrow(tempPopPop),
                                        pg$timestep), ]
  
  savePops <- tempPopPop$population
  tempPopPop <- select(tempPopPop, -population)
  
  
  #Calculate incidence per 100,000
  
  tempTotalInc <- data.frame(year = head(pg$years, -1), 
                         best = (1e5 * tempTotalInfections /
                           tempTotalPop)) %>%
    tbl_df()
  
  
  totalInc <- tempTotalInc %>%
    gather("sim", "incidence", 2:ncol(tempTotalInc)) %>%
    group_by(year) %>%
    summarise(min = min(incidence),
              max = max(incidence)) %>%
    ungroup() %>%
    mutate(best = tempTotalInc$best) %>%
    select(year, best, min, max)
  
  totalIncidencePlot <- indicatorPlot(totalInc, 
                                      ylabel = "Incidence per 100,000",
                                      range = FALSE, 
                                      xlimits = c(pg$start_year, 
                                                  pg$end_year, 20)) +
    ggtitle("Overall incidence")
  
  tempPopInc <- data.frame(year = rep(head(pg$years, -1), pg$npops),
                  population = savePops,
                    best = (1e5 * select(tempPopInfections, -population) /
                              tempPopPop)) %>%
    tbl_df()
  
  popIncidencePlot <- indicatorPlot(tempPopInc, 
                                    ylabel = "Incidence per 100,000",
                                    xlimits = c(pg$start_year, 
                                                pg$end_year, 20),
                                    range = FALSE, 
                                    groupPlot = "population") +
    ggtitle("Incidence by population")
  
  # Treatment 
  
  totalTreatments <- indicatorResults(pg, bestResults,
                                      "newTreatments",
                                      populations = "all", 
                                      annual = "sum")
  
  popTreatments <- indicatorResults(pg, bestResults, 
                                    "newTreatments", 
                                    annual = "sum")
  
  totalTreatmentPlot <- indicatorPlot(totalTreatments, 
                                  ylabel = "Number initiated treatment", 
                                  range = FALSE, 
                                  xlimits = c(pg$start_year, 
                                              pg$end_year, 20)) +
    ggtitle("Overall initiated treatment")
  
  popTreatmentPlot <- indicatorPlot(popTreatments, 
                                ylabel = "Number initiated treatment", 
                                xlimits = c(pg$start_year, 
                                            pg$end_year, 20),
                                range = FALSE, 
                                groupPlot = "population") +
    scale_colour_brewer(name = "Age Group", palette = "Set1",
                        labels = c("Age < 4", "Age 5 to 14", 
                                   "Age 15 to 44", "Age > 45")) +
    scale_fill_brewer(name = "Age Group", palette = "Set1",
                      labels = c("Age < 4", "Age 5 to 14", 
                                 "Age 15 to 44", "Age > 45")) +
    ggtitle("Initiated treatment by population")
  
  # HBV deaths 
  totalDeaths <- indicatorResults(pg, bestResults, 
                                  "newHBVdeaths", 
                                  populations = "all", 
                                  annual = "sum")
  
  popDeaths <- indicatorResults(pg, bestResults, 
                                "newHBVdeaths",
                                annual = "sum")
  
  totalDeathsPlot <- indicatorPlot(totalDeaths, 
                                   ylabel = "Number of deaths", 
                                   range = FALSE, 
                                   xlimits = c(pg$start_year, 
                                               pg$end_year, 20)) +
    ggtitle("Overall number of deaths")
  
  popDeathsPlot <- indicatorPlot(popDeaths, 
                                    ylabel = "Number of deaths", 
                                    xlimits = c(pg$start_year, 
                                                pg$end_year, 20),
                                    range = FALSE, 
                                    groupPlot = "population") +
    scale_colour_brewer(name = "Age Group", palette = "Set1",
                        labels = c("Age < 4", "Age 5 to 14", 
                                   "Age 15 to 44", "Age > 45")) +
    scale_fill_brewer(name = "Age Group", palette = "Set1",
                      labels = c("Age < 4", "Age 5 to 14", 
                                 "Age 15 to 44", "Age > 45")) +
    ggtitle("Deaths by population")
  
  # Put all plots in a grid
  print(plot_grid(totalPopPlot, 
                  agePopPlot, 
                  totalPrevPlot,
                  popPrevPlot, 
                  
                  totalNewInfectionsPlot,
                  totalCumInfectionsPlot,
                  popNewInfectionsPlot,
                  popCumInfectionsPlot,
                  
                  totalIncidencePlot,
                  popIncidencePlot, 
                  
                  totalTreatmentPlot, 
                  popTreatmentPlot, 
                  
                  totalDeathsPlot,
                  popDeathsPlot, 
                  ncol = 4, nrow = 4))
  
}

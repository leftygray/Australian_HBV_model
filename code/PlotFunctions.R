## Useful functions for HBV modelling and analysis

# R. T. Gray

# This script contains functions useful functions for plotting
# the HBV model results and for figure generation. 


# Extract annual results -------------------------------------------------

MidyearIndex <- function(n, timestep) {
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

SumYear <- function(vector, timestep) {
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
  yearFrame[, sumVar] <- SumYear(as.matrix(df[, sumVar]), timestep)
  yearFrame$year <- rep(head(years, -1), npops)
  
  if (!is.null(midVar)) {
    yearFrame[, midVar] <- df[MidyearIndex(nrow(df), timestep), midVar]
  }
  
  if (!is.null(divideVar)) {
    yearFrame[, divideVar] <- yearFrame[, sumVar] / 
      yearFrame[, midVar]
  }
  return(yearFrame)
}

# Organize results -------------------------------------------------------

FactorPop <- function(df, factorLabels) {
  df$population <- factor(df$population)
  levels(df$population) <- factorLabels
  
  return(df)
}

MidYearDf <- function(df, timstep, year, npops) {
  yearFrame <- df[seq(1, nrow(df), 1 / timestep), ]
  yearFrame <- df[MidyearIndex(nrow(df), timestep), ]
  yearFrame$year <- rep(head(years, -1), npops)
}

popResults <- function(pg, bestResults, paramResults,
                       populations = NULL, states = NULL,
                       range = FALSE) {
  # This function organizes and merges the population size results
  
  # First organize the best estimates
  bestPopSizes <- as.data.frame.table(bestResults$allPops)
  colnames(bestPopSizes) <- c("population", "state", "time_step", 
                              "popsize") 
  bestPopSizes$time_step <- as.numeric(bestPopSizes$time_step)
  
  bestPopSizes <- bestPopSizes %>% # add pts as well
    mutate(year = pg$start_year + (time_step - 1) * pg$timestep) %>%
    select(year, time_step, everything()) %>%
    select(-time_step) %>%
    arrange(population) %>%
    rename(best = popsize)
  
  # Now append each param set results
  popSizes <- bestPopSizes
  
  for (set in 1:pg$parameter_samples) {
    paramSet <- paste0("paramSet", toString(set))
    setPopSizes <- as.data.frame.table(paramResults[[paramSet]]$allPops)
    colnames(setPopSizes) <- c("population", "state", "time_step",
                               "popsize") 
    setPopSizes <- arrange(setPopSizes, population)
    popSizes[, paste0("set", toString(set))] <-  setPopSizes$popsize
    
  }

  # Extract the subresults we want
  if (!is.null(populations) || !is.null(states)) {
    if (length(populations) == 1 && populations == "all") {
      if (length(states) == 1 && states == "all") {
        popSizes <- popSizes %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      } else if (is.null(states)) {
        popSizes <- popSizes %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation, state) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      } else {
        popSizes <- popSizes %>%
          filter(state %in% states) %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      }
    } else if (is.null(populations)){
      if (length(states) == 1 && states == "all") {
        popSizes <- popSizes %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation, population) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      } else {
        popSizes <- popSizes %>%
          filter(state %in% states) %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation, population) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      }
    } else {
      if (length(states) == 1 && states == "all") {
        popSizes <- popSizes %>%
          filter(population %in% populations) %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      } else {
        popSizes <- popSizes %>%
          filter(state %in% states) %>%
          filter(population %in% populations) %>%
          gather("simulation", "popsize", 4:ncol(popSizes)) %>%
          group_by(year, simulation) %>%
          summarise(total_pop = sum(popsize)) %>%
          spread(simulation, total_pop)
      }
    }  
  }
  
  # Ungroup in preparation for next manipulation 
  popSizes <- ungroup(popSizes)
  bestValues <- popSizes$best
  
  if (range) {
    bestcol <- which(colnames(popSizes) == "best")
    
    popSizes <- gather(popSizes, "sim", "popsize", 
                       bestcol:ncol(popSizes)) 
    
    if (is.null(populations) && is.null(states)) {
      popSizes <- group_by(popSizes, year, population, state)
    } else if (is.null(populations)) {
      popSizes <- group_by(popSizes, year, population)
    } else if (is.null(states)) {
      popSizes <- group_by(popSizes, year, state)
    } else {
      popSizes <- group_by(popSizes, year)
    }
    
    popSizes <- popSizes %>%
      summarise(min = min(popsize),
                max = max(popsize)) %>%
      ungroup() 
    
    if (is.null(populations) && is.null(states)) {
      popSizes <- popSizes %>%
        arrange(population) %>%
        mutate(best = bestValues)
    } else {
      popSizes <- popSizes %>%
        mutate(best = bestValues)
    }
    
    if (is.null(populations) && is.null(states)) {
      popSizes <- select(popSizes, year, population, state,
                         best, everything())
    } else if (is.null(populations)) {
      popSizes <- select(popSizes, year, population, best, everything())
    } else if (is.null(states)) {
      popSizes <- select(popSizes, year, state, best, everything())
    } else {
      popSizes <- select(popSizes, year, best, everything())
    }
  }
  
  
  # Return the final population data frame
  return(tbl_df(popSizes))
}


indicatorResults<- function(pg, bestResults, paramResults,
                            indicator, populations = NULL, 
                            range = FALSE, 
                            annual = NULL) {
  # This function organizes and merges the specific
  # indicator results
  
  # annual = c("sum", "midyear", NULL)
  
  # First organize the best estimates
  bestEstimates <- tbl_df(as.data.frame(t(bestResults[[indicator]]))) %>%
    mutate(year = pg$pts) %>%
    select(year, everything()) %>%
    gather_("population", indicator, pg$population_names) %>%
    rename_("best" = indicator) 
  
  # Now append each param set results
  indicatorEstimates <- bestEstimates
  
  for (set in 1:pg$parameter_samples) {
    paramSet <- paste0("paramSet", toString(set))
    setIndicator <- 
      tbl_df(as.data.frame(t(paramResults[[paramSet]][[indicator]]))) %>%
      mutate(year = pg$pts) %>%
      select(year, everything()) %>%
      gather_("population", indicator, pg$population_names)
    
    # setIndicator <- arrange(setIndicator, population)
    indicatorEstimates[, paste0("set", toString(set))] <-  
      setIndicator[, indicator]
  }
  
  # Convert population into a factor and order the levels
  # indicatorEstimates$population <- factor(indicatorEstimates$population,
  #                                         levels = c("age0to4", "age5to14",
  #                                                    "age15to44", "age45"))
  
  # Extract the subresults we want
  if (length(populations) == 1 && populations == "all") {
    indicatorEstimates <- indicatorEstimates %>%
      gather("simulation", "estimate", 3:ncol(indicatorEstimates)) %>%
      group_by(year, simulation) %>%
      summarise(total = sum(estimate)) %>%
      spread(simulation, total)
  } else if (is.null(populations)){
    indicatorEstimates <- indicatorEstimates %>%
      gather("simulation", "estimate", 3:ncol(indicatorEstimates)) %>%
      group_by(year, simulation, population) %>%
      summarise(total = sum(estimate)) %>%
      spread(simulation, total)
  } else {
    indicatorEstimates <- indicatorEstimates %>%
      filter(population %in% populations) %>%
      gather("simulation", "estimate", 3:ncol(indicatorEstimates)) %>%
      group_by(year, simulation) %>%
      summarise(total = sum(estimate)) %>%
      spread(simulation, total)
  }  
  
  # Ungroup in preparation for next manipulation 
  indicatorEstimates <- ungroup(indicatorEstimates)
  if (is.null(populations)){
    indicatorEstimates <- arrange(indicatorEstimates, population)
  }
  
  # Extract best column location
  bestcol <- which(colnames(indicatorEstimates) == "best")
  
  # Apply annual manipulations
  if (!is.null(annual)) {
    if (annual == "sum") {
      tempEstimates <- indicatorEstimates[seq(1, nrow(indicatorEstimates), 
                                              1 / pg$timestep), ]
      tempEstimates[, bestcol:ncol(tempEstimates)] <- 
        apply(indicatorEstimates[, bestcol:ncol(indicatorEstimates)],
              2, function(x) SumYear(x, pg$timestep))
      
      if (is.null(populations)) {
        tempEstimates$year <- rep(head(pg$years, -1), pg$npops)
      } else {
        tempEstimates$year <- head(pg$years, -1)
      }
      
      indicatorEstimates <- tempEstimates
      
    } else if (annual == "midyear") {
      tempEstimates <- indicatorEstimates[seq(1, nrow(indicatorEstimates), 
                                              1 / pg$timestep), ]
      tempEstimates <- 
        indicatorEstimates[MidyearIndex(nrow(indicatorEstimates), 
                                        pg$timestep), ]
      if (is.null(populations)) {
        tempEstimates$year <- rep(head(pg$years, -1), pg$npops)
      } else {
        tempEstimates$year <- head(pg$years, -1)
      }
      
      indicatorEstimates <- tempEstimates
      
    }
  }
  
  # Store best values for later
  bestValues <- indicatorEstimates$best
  
  if (range) {
   
    indicatorEstimates <- gather(indicatorEstimates, "sim", "estimate", 
                       bestcol:ncol(indicatorEstimates)) 
    
    if (is.null(populations)) {
      indicatorEstimates <- indicatorEstimates %>%
        group_by(year, population) %>%
        summarise(min = min(estimate),
                  max = max(estimate)) %>%
        ungroup() %>%
        arrange(population) %>%
        mutate(best = bestValues) %>%
        select(year, population, 
               best, everything())
    } else {
      indicatorEstimates <- indicatorEstimates %>%
        group_by(year) %>%
        summarise(min = min(estimate),
                  max = max(estimate)) %>%
        ungroup() %>%
        mutate(best = bestValues) %>%
        select(year, best, everything())
    }
  }
  
  # Return the final population data frame
  return(indicatorEstimates)
  
}
  
  
# Plotting functions -----------------------------------------------------

indicatorPlot <- function(data, ylabel = NULL, range = TRUE,
                          xlimits = NULL, facetPlot = NULL, 
                          groupPlot = NULL) {

  # Set defaults
  if (is.null(ylabel)) {
    ylabel = indicator
  }
  
  if (is.null(xlimits)) {
    xlimits = c(min(data$year), max(data$year), 20)
  }
  
  if (is.null(facetPlot)) {
    plotFacets = FALSE
  } else {
    plotFacets = TRUE
  }
  
  if (is.null(groupPlot)) {
    plotGroups = FALSE
  } else {
    plotGroups = TRUE
  }
  
  # Create the plot
  if (range) {
    if (plotGroups) {
      plot <- ggplot(data = data, aes_string(x = "year", group = groupPlot)) +
        geom_ribbon(aes_string(ymin = "min", ymax = "max", 
                               fill = groupPlot), alpha = 0.4) +
        geom_line(aes_string(y = "best", colour = groupPlot)) +
        xlab("Year") + ylab(ylabel) +
        plotOpts + theme(legend.position = "right")
    } else {  
      plot <- ggplot(data = data, aes(x = year)) +
        geom_ribbon(aes(ymin = min, ymax = max),
                    fill = "blue", alpha = 0.4) +
        geom_line(aes(y = best), colour = "blue") +
        xlab("Year") + ylab(ylabel) +
        plotOpts
    }
  } else {
    if (plotGroups) {
      plot <- ggplot(data = data, aes_string(x = "year", group = groupPlot)) +
        geom_line(aes_string(y = "best", colour = groupPlot)) +
        xlab("Year") + ylab(ylabel) +
        plotOpts + theme(legend.position = "right")
    } else {  
      plot <- ggplot(data = data, aes(x = year)) +
        geom_line(aes(y = best), colour = "blue") +
        xlab("Year") + ylab(ylabel) +
        plotOpts
    }
  }
  
  plot <- plot + coord_cartesian(xlim = xlimits[1:2]) +
    scale_x_continuous(breaks = seq(xlimits[1], xlimits[2], 
                                    by = xlimits[3]))

  if (plotFacets) {
    plot <- plot + facet_wrap(as.formula(paste("~", facetPlot)),
                              ncol = 2, scales = "free") +
    theme(panel.margin = unit(2, "lines"))
  }
  
  # return the final plot
  return(plot)
}

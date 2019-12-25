

SIR <- function(initInfected, effectiveContacts, recoveryTime,
                               initSusceptible = 990,
                start = 0,
                finish = 30){
  
  library(deSolve)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  # Create the start time, finish time, and time step
  START  <- start
  FINISH <- finish
  STEP   <- 1 / 128
  
  # Create time vector
  simtime <- seq(START, FINISH, by = STEP)
  
  # Create stocks vector, with initial values
  stocks <- c(Susceptible         = initSusceptible,
              Infected            = initInfected,
              Recovered           = 0,
              Cumulative_Infected = initInfected)
  

  
  
  # Create auxiliaries vector, with values
  auxs    <- c(effectiveContacts = effectiveContacts, 
               recoveryTime = recoveryTime)
  
  model <- function (time, stocks, auxs) 
    with(as.list(c(stocks, auxs)), {
      
      population           <- Susceptible + Infected + Recovered
      probability          <- Susceptible / population
      contactsPerInfected  <- Infected * effectiveContacts
      #Flows
      RR                   <- Infected / recoveryTime
      IR                   <- contactsPerInfected  * probability
      # Net flows
      d_Susceptible_dt         <- -IR
      d_Infected_dt            <- IR - RR
      d_Recovered_dt           <- RR
      d_Cumulative_Infected_dt <- IR
      
      return(list(c(d_Susceptible_dt, d_Infected_dt, d_Recovered_dt,
                    d_Cumulative_Infected_dt), 
                  IR = IR, 
                  RR = RR, 
                  probability = probability, 
                  contactsPerInfected = contactsPerInfected, 
                  effectiveContacts = effectiveContacts, 
                  population = population, 
                  recoveryTime = recoveryTime))
    })
  
  # Call Solver, and store results in a data frame
  o <- data.frame(ode(y = stocks, times = simtime, func = model,
                      parms = auxs, method = "rk4"))
  
  tidy_sir <- pivot_longer(o, -time,
                            names_to = "variable", values_to = "value") %>%
    filter(variable %in% names(stocks))

  g <- ggplot(tidy_sir, aes(x = time, y = value)) +
    geom_line(aes(group = variable, colour = variable))
  
  incidence_data <- o %>% select(time, names(stocks)) %>% 
    filter(time - trunc(time) == 0) %>% 
    mutate(cumulativeInfected = Infected + Recovered,
           incidence = cumulativeInfected - lag(cumulativeInfected, 
                                                default = cumulativeInfected[1]),
           incidence = round(incidence)) %>% 
    select(time, incidence) %>% 
    filter(time != 0)
  
  list(incidence_data = incidence_data,
       model_df   = o,
       stocks_graph   = g)
}
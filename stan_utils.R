create_stan_file <- function(stan_text, filename) {
  file_connection <- file(filename)
  writeLines(stan_text, file_connection)
  close(file_connection)
}

# n_constants is the number of unknown parameters (except init values) in the
# ODE Model

fit_stan_model <- function(n_constants, incidence_data,
                           warm_up, iterations, seed = NULL,
                           adapt_delta = 0.8) {
  
  stan_d = list(n_obs    = length(incidence_data),
                n_params = n_constants, 
                n_difeq  = 3, # number of differential equations
                y        = incidence_data,
                t0       = 0,
                ts       = 1:30)
  
  # Test / debug the model:
  test <- stan("SIR.stan",
               data = stan_d,
               chains = 1, 
               iter = 10,
               verbose = FALSE,
               refresh = 0)
  
  if(is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  
  # Fit and sample from the posterior
  stan_fit <- stan(fit    = test,
                   data   = stan_d,
                   chains = 3,
                   warmup = warm_up,
                   iter   = iterations,
                   cores  = 3,
                   seed    = seed,
                   control = list(adapt_delta = adapt_delta))
}

generate_stan_model_three_params <- function(prior_susceptible, 
                                             likelihood = "normal") {
  
  if(likelihood == "normal") {
    
    model <- paste("
model {
  params[1] ~ uniform(0, 500);
  params[2] ~ uniform(0, 20);",
                   prior_susceptible,
                   "sigma ~ cauchy( 0 , 2 );
  y ~ normal(incidence, sigma); 
}
", sep = "\n  ")
  }

  if(likelihood == "poisson") {
    
    model <- paste("
model {
  params[1] ~ uniform(0, 500);
  params[2] ~ uniform(0, 20);",
                   prior_susceptible,
  "y ~ poisson(incidence);
}
  ", sep = "\n  ")
  }

  model
}



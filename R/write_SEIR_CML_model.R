write_SEIR_CML_model <- function(matrix_type, filename, stocks, params_prior) {
  source("./stan_utils.R")
  
  if(matrix_type == "A") {
    model_file   <- "./deterministic_models/cml_4_cohorts_SEIR_matrix_A.stmx"
    ODE_fun_name <- "SEIR_cml_matrix_A"
  }
  
  if(matrix_type == "B") {
    model_file   <- "./deterministic_models/cml_4_cohorts_SEIR_matrix_B.stmx"
    ODE_fun_name <- "SEIR_cml_matrix_B"
  }
  
  mdl            <- read_xmile(model_file)
  constant_names <- map_chr(mdl$description$constants, "name")
  # serial interval parameters
  si_pars        <- c("recovery_time", "latent_period")
  params         <- constant_names[!constant_names%in% si_pars ] %>% sort()
  
  override_consts <- list(
    list(name = "recovery_time", value = 2),
    list(name = "latent_period", value = 1))
  
  stan_fun   <- create_stan_function(model_file, ODE_fun_name, 
                                     pars = params,
                                     override.consts = override_consts)
  
  stan_data  <- paste(
    "data {",
    "  int<lower = 1> n_obs; // Number of weeks sampled",
    "  int<lower = 1> n_params; // Number of model parameters",
    "  int<lower = 1> n_difeq; // Number of differential equations in the system",
    "  real y1[n_obs];",
    "  real y2[n_obs];",
    "  real y3[n_obs];",
    "  real y4[n_obs];",
    "  real t0; // Initial time point (zero)",
    "  real ts[n_obs]; // Time points that were sampled",
    "}", sep = "\n")
  
  stan_td     <- generate_transformed_data_block()
  
  stan_params <- paste(
    "parameters {",
    "  real<lower = 0> params[n_params]; // Model parameters",
    "  real<lower = 0> sigma1;",
    "  real<lower = 0> sigma2;",
    "  real<lower = 0> sigma3;",
    "  real<lower = 0> sigma4;",
    "}", sep = "\n")
  
  counter <- 0
  
  stan_init_stocks <- sapply(stocks, function(stock) {
    counter <<- counter + 1
    paste0("  y0[", counter, "] = ", stock, ";")
  }) %>% paste(collapse = "\n")
  
  fun_exe_line <- paste0("  y_hat = integrate_ode_rk45(",
                         ODE_fun_name,", y0, t0, ts, params, x_r, x_i);")
  
  stan_tp_cml <- paste(
    "transformed parameters{",
    "  real y_hat[n_obs, n_difeq]; // Output from the ODE solver",
    "  real y0[n_difeq]; // Initial conditions",
    "  real cum_inc1[n_obs];",
    "  real cum_inc2[n_obs];",
    "  real cum_inc3[n_obs];",
    "  real cum_inc4[n_obs];",
    stan_init_stocks,
    fun_exe_line,
    "  for (i in 1:n_obs) {",
    "    cum_inc1[i] = y_hat[i, 17];",
    "    cum_inc2[i] = y_hat[i, 18];",
    "    cum_inc3[i] = y_hat[i, 19];",
    "    cum_inc4[i] = y_hat[i, 20];",
    "  }",
    "}", sep = "\n")
  
  stan_model_text <- c(
    params_prior,
    "  sigma1 ~ cauchy( 0 , 2 )",
    "  sigma2 ~ cauchy( 0 , 2 )",
    "  sigma3 ~ cauchy( 0 , 2 )",
    "  sigma4 ~ cauchy( 0 , 2 )",
    "  y1     ~ normal(cum_inc1, sigma1)",
    "  y2     ~ normal(cum_inc2, sigma2)",
    "  y3     ~ normal(cum_inc3, sigma3)",
    "  y4     ~ normal(cum_inc4, sigma4)"
  )
  stan_mdl  <- generate_model_text(stan_model_text)
  
  stan_text   <- paste(stan_fun, stan_data, stan_td, stan_params,
                       stan_tp_cml, stan_mdl, sep = "\n")
  
  create_stan_file(stan_text, filename)
  
  list(params = params)
}
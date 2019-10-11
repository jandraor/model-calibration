source("./deterministic_models/SIR.R")

library(dplyr)
library(purrr)

real_data <- incidence_data_SIR(10, 1, 2, 990)

fits <- c("S0 ~ normal(990, 1);",
          "S0 ~ normal(791, 10);",
          "S0 ~ normal(791, 250);")

fitted_df <- map2_df(fits, stan_fits, function(fit, stan_fit) {
  stan_df <- as.data.frame(stan_fit)
  ec      <- mean(stan_df$`params[1]`)
  print(ec)
  rt      <- mean(stan_df$`params[2]`)
  print(rt)
  S0      <- mean(stan_df$S0)
  print(S0)
  cat("\n R_nought:", ec * rt, "\n")
  fit_data <- incidence_data_SIR(10, ec, rt, S0) %>% 
    mutate(type = fit)
})



dataset <- bind_rows(real_data, fitted_df)

ggplot(fitted_df, aes(x = time, y = incidence)) +
  geom_line(aes(group = type), alpha = 0.5) +
  facet_wrap(~type) +
  geom_line(data = real_data)
  
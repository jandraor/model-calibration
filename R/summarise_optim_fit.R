summarise_optim_fit <- function(fit, incidence_data, conceptual_matrix, pop_cohorts) {
  translation_df <- tibble(
    index_group = 1:4,
    var_incidence = paste0("incidence", 1:4),
    age_group = age_groups <- c("00-04", "05-14", "15-44", "45+"))
  
  actual_df <- incidence_data %>% left_join(translation_df) %>% 
    select(-var_incidence, -index_group) %>% mutate(source = "syn_data")
  
  comparison_df <- bind_rows(actual_df, fit$fit_df)
  
  g_comparison <- ggplot(comparison_df, aes(x = time, y = incidence)) +
    geom_line(aes(group = source, colour = source)) +
    scale_colour_manual(values = c("steelblue", "grey")) +
    facet_wrap(~age_group) +
    theme_test()
  
  WAIFW_estimates <- sapply(
    conceptual_matrix, function(Bij, row) row[[Bij]], row = optim_fit$fit$par) %>%
    matrix(nrow = 4, byrow = T)
  
  age_groups <- c("Age.0.4", "Age.5.14", "Age.15.44", "Age.45.over")
  colnames(WAIFW_estimates) <- rownames(WAIFW_estimates) <- age_groups
  g_WAIFW <- draw_WAIFW(WAIFW_estimates, "title")
  
  next_gen_matrix <- WAIFW_estimates / 1e5 * pop_cohorts * 2
  eigensystem     <- eigen(next_gen_matrix)
  R_nought        <- max(abs(eigensystem$values))
  
  list(g_comparison = g_comparison,
       g_WAIFW = g_WAIFW,
       R_nought = R_nought)
}
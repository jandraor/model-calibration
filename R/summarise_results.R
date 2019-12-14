summarise_results <- function(fit, conceptual_matrix, incidence_df, pop_sizes, specs_list, unique_params) {
  source("./graphs.R")
  
  posterior_df <- as.data.frame(fit)
  
  mean_incidences <- posterior_df %>% select(contains("incidence")) %>% 
    colMeans()
  
  names_mi <- names(mean_incidences)
  
  age_groups <- c("Age.0.4", "Age.5.14", "Age.15.44", "Age.45.over")
  
  sim_data <- map_df(age_groups, function(ag) {
    i <- which(ag == age_groups)
    regex <- stringr::regex(paste0("incidence", i, "\\["))
    data.frame(week = 1:50, value = mean_incidences[grepl(regex, names_mi)], 
               variable = ag)
  }) %>% mutate(source = "sim data")
  
  samples                 <- rstan::extract(fit)
  param_samples           <- samples$params
  colnames(param_samples) <- unique_params
  param_medians           <- apply(param_samples, 2, function(col) median(col))
  
  WAIFW_medians <- sapply(
    conceptual_matrix, function(Bij, row) row[[Bij]], row = param_medians) %>%
    matrix(nrow = 4, byrow = T)
  
  normalised_WAIFW           <- WAIFW_medians
  colnames(normalised_WAIFW) <- rownames(normalised_WAIFW) <- age_groups
  
  g_WAIFW <- draw_WAIFW(normalised_WAIFW, "title")
  
  #============================================================================
  
  mean_incidences <- posterior_df %>% select(contains("incidence")) %>% 
    colMeans()
  
  names_mi <- names(mean_incidences)
  
  age_groups <- c("Age.0.4", "Age.5.14", "Age.15.44", "Age.45.over")
  
  translation_df <- data.frame(index_group = 1:4, variable = age_groups)
  
  sim_data <- map_df(age_groups, function(ag) {
    i <- which(ag == age_groups)
    regex <- stringr::regex(paste0("incidence", i, "\\["))
    data.frame(time = 1:50, value = mean_incidences[grepl(regex, names_mi)], 
               variable = ag)
  }) %>% mutate(source = "sim data")
  
  real_data <- incidence_df %>% rename(value = incidence) %>% 
    mutate(source = "syn data") %>% left_join(translation_df, by = "index_group")
  
  comparison_data <- bind_rows(sim_data, real_data)
  
  g_comparison <- ggplot(comparison_data, aes(x = time, y = value)) +
    geom_line(aes(group = source, colour = source)) +
    scale_colour_manual(values = c("lightgrey", "blue")) +
    facet_wrap(~variable) +
    theme_test()
  
  #============================================================================
 
  r_noughts <- sapply(1:nrow(param_samples), function(i) {
    row <- param_samples[i, ]
    WAIFW <- sapply(conceptual_matrix, 
                    function(Bij, row) row[[Bij]], row = row) %>%
      matrix(nrow = 4, byrow = TRUE)
    next_gen_matrix <- WAIFW / 1e5 * pop_sizes * 2
    eigensystem     <- eigen(next_gen_matrix)
    max(abs(eigensystem$values))
  })
  
  g_rNougths <- draw_density(r_noughts, specs_list)
  
  list(g_WAIFW = g_WAIFW,
       g_comparison = g_comparison,
       g_rNougths = g_rNougths,
       mean_rNought = mean(r_noughts))
}
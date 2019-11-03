draw_multiple_densities <- function(densities_df, stats_df,
                                    x_pos, text_size = 5) {
  g <- ggplot(densities_df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
    facet_wrap(~ param, scales = "free") +
    scale_fill_manual(values = c(NA, "lightgrey"))  +
    geom_vline(aes(xintercept = mean_value), stats_df, color = "#1261A0",
               linetype ="dashed", size = 1) +
    geom_text(aes(x = x_pos, y = y_pos_mean, label = mean_label),
              stats_df, colour = "#1261A0", size = text_size) +
    geom_text(aes(x = x_pos, y = y_pos_median, label = median_label),
              stats_df, colour = "#1261A0", size = text_size) +
    geom_text(aes(x = x_pos, y = y_pos_interval, label = interval_label),
               stats_df, colour = "grey", size = text_size) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

generate_graph_inputs <- function(stan_fit, x_pos, y_pos, y_pos_median, 
                                  y_pos_interval, pars, rename_pars = NULL) {
  
  posterior_df <- as.data.frame(stan_fit) 
  
  for(rename_list in rename_pars) {
    posterior_df <- posterior_df %>% 
      rename(!!rename_list$new := !!rename_list$old)
  }
  
  if("recoveryTime" %in% pars) {
    posterior_df <- mutate(posterior_df, recoveryTime = 1 / recoveryProportion)
  }
  
  if("latent_period" %in% pars) {
    posterior_df <- mutate(posterior_df, 
                           latent_period = 1 / incubation_proportion)
  }
     
  params_df <- select(posterior_df, pars)
  
  credible_intervals <- apply(params_df, 2, HPDI, prob = 0.95) %>% t() %>% 
    as.data.frame() %>% 
    rename(lower_interval = "|0.95", upper_interval = "0.95|") %>% 
    mutate(param = rownames(.))
  
  means   <- apply(params_df, 2, mean)
  medians <- apply(params_df, 2, median)
  
  stats_df <- data.frame(stringsAsFactors = FALSE,
                         param = names(means), 
                         mean_value = means,
                         median_value = medians,
                         lower_interval = credible_intervals$lower_interval,
                         upper_interval = credible_intervals$upper_interval) %>% 
    mutate(mean_label     = paste0("mean   = ", round(mean_value, 3)),
           median_label   = paste0("median = ", round(median_value, 3)),
           interval_label = paste0("[ ", round(lower_interval, 3), ", ", 
                                   round(upper_interval, 3), " ]"),
           x_pos          = x_pos,
           y_pos_mean     = y_pos,
           y_pos_median   = y_pos_median,
           y_pos_interval = y_pos_interval)
  
  densities <- apply(params_df, 2, density) %>% lapply(function(densityObj){
    data.frame(x = densityObj$x, y = densityObj$y)
  }) %>% bind_rows(.id = "param") %>% 
    inner_join(credible_intervals) %>% 
    mutate(area = x >= lower_interval & x <= upper_interval)
  
  list(densities = densities,
       stats_df  = stats_df)
}
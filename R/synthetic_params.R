produce_synthetic_params <- function() {
  
  # Data obtained S5 table from Mossong et al (2018)
  file                     <- "./data/contact_matrix_Finland.csv"
  raw_matrix               <- read_csv(file)
  raw_M_matrix             <- as.matrix(raw_matrix[, -1]) %>% t()
  age_cohorts              <- raw_matrix[, 1] %>% pull() 
  colnames(raw_M_matrix)   <- age_cohorts   
  row.names(raw_M_matrix)  <- age_cohorts
  
  g_raw_M_matrix <- draw_WAIFW(raw_M_matrix, "Finland's social contact matrix")
  
  raw_data <- read_csv("./data/Finland_population_data.csv") 
  
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  
  fin_pop <- raw_data %>% 
    select(AgeGrp, AgeGrpStart, PopTotal) %>% 
    mutate(AgeGrp = ifelse(AgeGrp == "0-4",  "00-04", AgeGrp),
           AgeGrp = ifelse(AgeGrp == "5-9",  "05-09", AgeGrp),
           AgeGrp = ifelse(AgeGrpStart >= 70, "70+", AgeGrp),
           AgeGrpStart = ifelse(AgeGrpStart >= 70, "70", AgeGrpStart)) %>% 
    group_by(AgeGrp, AgeGrpStart) %>% summarise(PopTotal = sum(PopTotal)) %>% 
    ungroup() %>% arrange(as.numeric(AgeGrpStart)) %>% 
    select(AgeGrp, PopTotal) 
  
  grouping_table <- data.frame(stringsAsFactors = FALSE,
                               five_yr_cohort = unique(fin_pop$AgeGrp),
                               group = c(age_groups[1], 
                                         rep(age_groups[2], 2), 
                                         rep(age_groups[3], 6),
                                         rep(age_groups[4], 6)))
  
  fin_pop <- fin_pop %>% 
    left_join(grouping_table, by = c("AgeGrp" = "five_yr_cohort"))
  
  # Aggregated population by defined age groups
  agg_pop <- fin_pop %>% group_by(group) %>% 
    summarise(population = sum(PopTotal)) %>% ungroup()
  
  total_pop <- sum(agg_pop$population)
  
  agg_pop$proportion <- agg_pop$population / total_pop
  
  agg_pop$trans_factor <- with(agg_pop, 1 / proportion)
  
  contacts_df <- melt(raw_M_matrix) %>% 
    rename(contactor = Var1, contactee = Var2) %>% 
    select(contactor, contactee, value) %>% 
    left_join(fin_pop, by = c("contactor" = "AgeGrp")) %>% 
    rename(contactor_group = group)
  
  grouped_contacts_df <- contacts_df %>% 
    group_by(contactor_group, contactee) %>% 
    mutate(subpopulation = sum(PopTotal),
           weight = PopTotal / subpopulation,
           contacts_weighted = value * weight) %>% 
    summarise(contacts = sum(contacts_weighted)) %>% 
    left_join(grouping_table, by = c("contactee" = "five_yr_cohort")) %>% 
    rename(contactee_group = group) %>% 
    group_by(contactor_group, contactee_group) %>% 
    summarise(contacts = sum(contacts))
  
  subtitle <- "Aggregated contacts"
  
  g_M_matrix <- ggplot(
    data = grouped_contacts_df, aes(
      x = contactor_group,
      y = ordered(contactee_group, levels = rev(sort(unique(contactee_group)))),
      fill = contacts)) + 
    geom_tile() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    geom_text(aes(label = round(contacts)), colour = "white", size = 2) +
    theme_minimal() + 
    labs(y ="", x = "",subtitle = subtitle) +
    theme(legend.position = "none",
          plot.subtitle = element_text(color = "#404040", size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
  
  wider_df <- pivot_wider(grouped_contacts_df, names_from = contactor_group, 
                          values_from = contacts)
  # Social contact matrix
  M_matrix <- wider_df[, -1] %>% as.matrix()
  rownames(M_matrix) <- colnames(M_matrix)
  
  C_matrix <- M_matrix * agg_pop$trans_factor
  
  # Symmetric normalised age-specific contact rate
  symmetric_C_matrix <- (C_matrix + t(C_matrix)) / 2
  
  # Social contact matrix corrected for reciprocity
  corrected_M_matrix <- symmetric_C_matrix * (1 / agg_pop$trans_factor) 
  
  infectious_period <- 2 # days  
  infectivity       <- 0.1
  syn_pop_size      <- 1e4
  
  # Normalised age-specific effective contact rate
  necr_matrix <- symmetric_C_matrix * infectivity
  
  synthetic_population <- agg_pop %>% 
    mutate(total_pop = sum(population),
           proportion = population / total_pop,
           syn_pop = round(syn_pop_size * proportion)) %>% 
    select(group, syn_pop)
  
  #=============================================================================
  # Theoretical R0
  #=============================================================================
  next_generation_matrix <- corrected_M_matrix * infectious_period * infectivity
  eigensystem            <- eigen(next_generation_matrix)
  theoretical_R0         <-  max(Re(eigensystem$values))
  
  syn_WAIFW <- necr_matrix / syn_pop_size 
  
  g_syn_WAIFW <- draw_WAIFW(syn_WAIFW * 1e5, "")
  
  list(g_original_contacts = g_raw_M_matrix,
       g_M_matrix = g_M_matrix,
       synthetic_WAIFW = syn_WAIFW,
       g_syn_WAIFW = g_syn_WAIFW,
       syn_pop = synthetic_population,
       theoretical_R0 = theoretical_R0)
}


  
  
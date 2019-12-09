produce_synthetic_params <- function() {
  library(readr)
  source("./graphs.R")
  
  # Data obtained S5 table from Mossong et al (2018)
  file                     <- "./data/contact_matrix_Finland.csv"
  raw_matrix               <- read_csv(file)
  contact_matrix_Fin       <- as.matrix(raw_matrix[, -1]) %>% t()
  age_cohorts              <- raw_matrix[, 1] %>% pull() 
  colnames(contact_matrix_Fin)  <- age_cohorts   
  row.names(contact_matrix_Fin) <- age_cohorts
  
  g_contacts_Fin <- draw_WAIFW(contact_matrix_Fin, "Finland's contact matrix")
  
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
  
  
  
  contacts_df <- melt(contact_matrix_Fin) %>% rename(from = Var2, to = Var1) %>% 
    select(from, to, value) %>% 
    left_join(fin_pop, by = c("from" = "AgeGrp")) %>% 
    rename(from_group = group)
  
  grouped_contacts_df <- contacts_df %>% group_by(from_group, to) %>% 
    mutate(subpopulation = sum(PopTotal),
           weight = PopTotal / subpopulation,
           contacts_weighted = value * weight) %>% 
    summarise(contacts = sum(contacts_weighted)) %>% 
    left_join(grouping_table, by = c("to" = "five_yr_cohort")) %>% 
    rename(to_group = group) %>% group_by(from_group, to_group) %>% 
    summarise(contacts = sum(contacts))
  
  subtitle <- "Aggregated contacts"
  
  g_grouped_contacts <- ggplot(
    data = grouped_contacts_df, aes(
      x = from_group,
      y = ordered(to_group, levels = rev(sort(unique(to_group)))),
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
  
  wider_df <- pivot_wider(grouped_contacts_df, names_from = from_group, values_from = contacts)
  grouped_contact_matrix            <- wider_df[ , -1] %>% as.matrix()
  rownames(grouped_contact_matrix ) <- colnames(grouped_contact_matrix)
  
  #effective group contact matrix
  infectivity <- 0.1
  egc_matrix  <- grouped_contact_matrix * infectivity
  
  synthetic_population <- fin_pop %>% group_by(group) %>% 
    summarise(population = sum(PopTotal)) %>% ungroup() %>% 
    mutate(total_pop = sum(population),
           proportion = population / total_pop,
           syn_pop = round(10000 * proportion)) %>% 
    select(group, syn_pop)
  
  eigensystem    <- eigen(egc_matrix * 2)
  theoretical_R0 <-  max(abs(eigensystem$values))
  
  syn_WAIFW <- egc_matrix / synthetic_population$syn_pop
  
  g_syn_WAIFW <- draw_WAIFW(syn_WAIFW * 1e5, "")
  
  list(g_original_contacts = g_contacts_Fin,
       g_aggregated_contacts = g_grouped_contacts,
       synthetic_WAIFW = syn_WAIFW,
       g_syn_WAIFW = g_syn_WAIFW,
       syn_pop = synthetic_population,
       theoretical_R0 = theoretical_R0)
}


  
  
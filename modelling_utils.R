generate_infection_rate <- function(index_S, indexes_I) {
  library(magrittr)
  sapply(indexes_I, function(index_I){
    paste0("B", index_S, index_I, "*", "I", index_I)
  }) %>% paste(collapse = "+")
}
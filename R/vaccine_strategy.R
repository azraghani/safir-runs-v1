get_vaccine_strategy <- function(days_to_vacc_start, doses_per_day, time_period, max_coverage, age_groups_covered, age_groups_covered_d3, vaccine_doses, pop, vacc_per_week, t_d3, t_10y_start){
  
    # we want to simulate a country having vaccinated all of the eligible population with the first and second dose before switching to a booster dose
    vaccine_set <- c(rep(0, days_to_vacc_start), rep(doses_per_day, time_period - days_to_vacc_start))
    
    # if vaccinating children < 10 years, want to wait until 1 Jan 2022
    if (age_groups_covered >= 16){
      # how many days to vaccinate population >10 years with 2 doses?
      population_10y <- sum(pop$n[3:17]) / sum(pop$n)
      days_vacc_10y <- floor(population_10y / (vacc_per_week/7) * max_coverage * 2)
     vaccine_set <- c(rep(0, days_to_vacc_start),
                       rep(doses_per_day, days_vacc_10y),
                       rep(0, max(t_10y_start - days_to_vacc_start - days_vacc_10y,0)),
                       rep(doses_per_day, max(time_period - t_10y_start,0)))
    }
    
    
    vaccine_coverage_strategy <- list()
    
    if (vaccine_doses == 2) {
      vaccine_coverage_strategy[[1]] <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = max_coverage)[1:age_groups_covered,]
      vaccine_coverage_strategy[[2]] <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = max_coverage)[1:age_groups_covered,]
    } else if (vaccine_doses == 3) {
      vaccine_coverage_strategy[[1]] <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = max_coverage)[1:age_groups_covered,]
      vaccine_coverage_strategy[[2]] <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = max_coverage)[1:age_groups_covered,]
      vaccine_coverage_strategy[[3]] <- nimue::strategy_matrix(strategy = "Elderly",max_coverage = max_coverage)[1:age_groups_covered_d3,]}
    
    next_dose_priority <- matrix(data = 1, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))
    
    # don't prioritise anyone for a booster until all population has received primary series
    if (vaccine_doses == 3){
      next_dose_priority[2,] <- 0
    }
  
  return(list(vaccine_set = vaccine_set, vaccine_coverage_strategy = vaccine_coverage_strategy, next_dose_priority = next_dose_priority, t_d3 = t_d3))
  
}

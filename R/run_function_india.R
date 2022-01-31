# main running function
run_scenario <- 
  function(fit, 
           scenario=1, 
           target_pop,
           hs_constraints,
           vacc_start, 
           vaccine_doses,
           vaccine,
           max_coverage,
           age_groups_covered,
           dt = 0.25,
           repetition = 1,
           seeding_cases = 10,
           vfr,
           vfr_time1,
           vfr_time2,
           variant_fold_reduction = 1,
           dose_3_fold_increase,
           age_groups_covered_d3,
           vacc_per_week,
           name = "scenario1",
           ab_model_infection = TRUE,
           std10 = 0.44,
           std10_infection = 0.44,
           t_d3,
           mu_ab_d1,
           mu_ab_d2,
           k,
           hl_s,
           hl_l,
           period_s,
           period_l,
           ab_50,
           ab_50_severe,
           mu_ab_infection,
           max_Rt,
           max_Rt_omicron,
           hosp_scal_omicron,
           ICU_scal_omicron){
    
    # set up transmission
    
    # start of simulation
    #R0_t0 <- as.Date(x = "2/1/2020", format = "%m/%d/%Y") [JT changed to below for excess fits]
    R0_t0 <- as.Date(x = Rt_overtime$date[1], format = "%m/%d/%Y") #take first date of excess fits 2019-12-14
    
    # three pulses to time waves correctly 
    R0_t1 <- as.Date(x = "3/1/2020", format = "%m/%d/%Y")
    R0_t2 <- as.Date(x = "11/1/2020", format = "%m/%d/%Y")
    R0_t3 <- as.Date(x = "6/1/2021", format = "%m/%d/%Y")
    
    # omicron introduction
    R0_t4 <- as.Date(x = vfr_time1, format = "%m/%d/%Y")
    
    # end of simulation
    tmax_date <- as.Date(x = "01/01/2022", format = "%m/%d/%Y")
    
    time_period <- as.integer(difftime(tmax_date, R0_t0 - 1))
 
    # get index for 1 Jan 2022, as don't want to vaccinate any children <10y before this time
    t_10y_start <- as.integer(difftime(as.Date("01/01/2022", format = "%m/%d/%Y"), R0_t0-1))
 
    # read in saved R profile [JT- replaced with below]
    # R_profile <- readRDS("data/R_profile_india.RDS")
    # R_profile$date <- as.Date(R_profile$date, format = "%d/%m/%Y")
    # add_R_profile <- data.frame(date=c(as.Date(x = vfr_time1, format = "%m/%d/%Y"), as.Date(x = vfr_time2, format = "%m/%d/%Y")), 
    #                             Rt=c(max_Rt,max_Rt_omicron))
    # R_profile <- rbind(R_profile,add_R_profile)
    # 
    # dates <- R_profile$date
    # rt <- R_profile$Rt
    # 
    # rt_out <- safir::interpolate_rt(dates = dates, rt = rt, max_date = tmax_date) 
    
    ####### use R profile from excess fits
    # read in saved R profile
    R_profile <- Rt_overtime
    R_profile$date <- as.Date(Rt_overtime$date, format = "%d/%m/%Y")
    missing_dates <- seq(as.Date(vfr_time1, format = "%m/%d/%Y"),as.Date(vfr_time2, format = "%m/%d/%Y"),by="days")
    missing_dates <- missing_dates[2:(length(missing_dates)-1)]
    missing_rt <- rep(NA, length(missing_dates))
    add_R_profile <- data.frame(date=c(as.Date(x = vfr_time1, format = "%m/%d/%Y"), missing_dates, as.Date(x = vfr_time2, format = "%m/%d/%Y"),
                                       as.Date(x = vfr_time2, format = "%m/%d/%Y")+1), 
                                Rt=c(max_Rt,missing_rt, max_Rt_omicron, max_Rt_omicron))
    
    R_profile <- rbind(R_profile,add_R_profile)
    
    dates <- R_profile$date
    rt <- R_profile$Rt
    rt <- na.approx(R_profile$Rt) #linear interpolation to fill in NAs
    
    rt_out<- list("Rt" = rt, "Rt_tt" = 1:length(rt))
    ###########################
    
    # daily per-capita prob of external infection
    lambda_external <- rep(0.0000001, time_period)
    
    # first pulse, spread out hazard of 0.001 over 10 days right before 1st wave
    t_spread <- 10
    lambda_tt <- as.integer(difftime(R0_t1, R0_t0 - 1))
    lambda_tt <- seq(from = lambda_tt - t_spread/2, to = lambda_tt + t_spread/2, by = 1)
    lambda_external[lambda_tt] <- dnorm(x = lambda_tt, mean = lambda_tt[t_spread/2+1], sd = 3) * 0.001
       
    # first pulse, spread out hazard of 0.001 over 10 days right before 2nd wave
    t_spread <- 10
    lambda_tt <- as.integer(difftime(R0_t2, R0_t0 - 1))
    lambda_tt <- seq(from = lambda_tt - t_spread/2, to = lambda_tt + t_spread/2, by = 1)
    lambda_external[lambda_tt] <- dnorm(x = lambda_tt, mean = lambda_tt[t_spread/2+1], sd = 3)
    lambda_external[lambda_tt] <- lambda_external[lambda_tt] / sum(lambda_external[lambda_tt]) * 0.001
       
    # second pulse, spread out hazard of 0.001 over 10 days right before 2nd wave
    t_spread <- 10
    lambda_tt <- as.integer(difftime(R0_t3, R0_t0 - 1))
    lambda_tt <- seq(from = lambda_tt - t_spread/2, to = lambda_tt + t_spread/2, by = 1)
    lambda_external[lambda_tt] <- dnorm(x = lambda_tt, mean = lambda_tt[t_spread/2+1], sd = 3)
    lambda_external[lambda_tt] <- lambda_external[lambda_tt] / sum(lambda_external[lambda_tt]) * 0.001
       
       
    # third pulse, spread out hazard of 0.001 over 10 days right before omicron wave
    t_spread <- 10
    lambda_tt <- as.integer(difftime(R0_t4, R0_t0 - 1))
    lambda_tt <- seq(from = lambda_tt - t_spread/2, to = lambda_tt + t_spread/2, by = 1)
    lambda_external[lambda_tt] <- dnorm(x = lambda_tt, mean = lambda_tt[t_spread/2+1], sd = 3)
    lambda_external[lambda_tt] <- lambda_external[lambda_tt] / sum(lambda_external[lambda_tt]) * 0.001
    
    # timing of vaccine introduction
    
    vacc_start <- as.Date(x = vacc_start, format = "%m/%d/%Y")
    days_to_vacc_start <- as.integer(difftime(vacc_start, R0_t0))
    
    # BASE PARAMETERS
    
    # Population and mixing
    rep_country <- "India"
    iso3c <- countrycode(rep_country, origin = "country.name", destination = "iso3c")
    pop <- squire::get_population(country = rep_country)
    pop_standardise <- target_pop / sum(pop$n)
    pop$n <- as.integer(pop$n * pop_standardise)
    contact_mat <- squire::get_mixing_matrix(country = rep_country)
    
    # Hospital capacity
    hc <- get_capacity(country = rep_country, pop = pop$n, hs_constraints = hs_constraints)
    
    # Poorer health outcomes for LMICs
    pnsdt <- squire:::probs$prob_non_severe_death_treatment
    
    # base parameters
    parameters <- safir::get_parameters(
      population = pop$n,
      contact_matrix_set = contact_mat,
      iso3c = iso3c,
      R0 = rt_out$Rt,
      tt_R0 = rt_out$Rt_tt,
      time_period = time_period,
      dt = dt,
      hosp_bed_capacity = hc$hosp_beds,
      ICU_bed_capacity = hc$ICU_beds,
      prob_non_severe_death_treatment = pnsdt,
      seeding_cases = seeding_cases,
      lambda_external = lambda_external
    )
    
    # --------------------------------------------------------------------------------
    # get vaccine parameters
    # --------------------------------------------------------------------------------
    # doses available each day
    doses_per_day <- floor(sum(pop$n) * vacc_per_week / 7)
    
    vaccine_out <- get_vaccine_strategy(days_to_vacc_start, doses_per_day, time_period, max_coverage, age_groups_covered, 
                                        age_groups_covered_d3, vaccine_doses, pop, vacc_per_week, t_d3, t_10y_start)
    
    vaccine_set <- vaccine_out$vaccine_set
    vaccine_coverage_strategy <- vaccine_out$vaccine_coverage_strategy
    next_dose_priority <- vaccine_out$next_dose_priority
    t_d3 <- vaccine_out$t_d3
    
    # profiles and dosing
    ## note that we later modify the decay rate vector that has been calculated here
    
    vax_pars <- get_vaccine_pars(vaccine = vaccine,
                                 mu_ab_d1 = mu_ab_d1,
                                 mu_ab_d2 = mu_ab_d2,
                                 k = k,
                                 hl_s = hl_s,
                                 hl_l = hl_l,
                                 period_s = period_s,
                                 period_l = period_l,
                                 dose_3_fold_increase = dose_3_fold_increase,
                                 ab_50 = ab_50,
                                 ab_50_severe = ab_50_severe,
                                 vaccine_doses = vaccine_doses,
                                 std10 = std10,
                                 t_d3 = t_d3)
    
    # dosing
    if (vaccine_doses == 2) {dose_period <- c(NaN, 28)}
    if (vaccine_doses == 3) {dose_period <- c(NaN, 28, (t_d3 + 28))}
    
    # combine parameters and verify
    parameters <- make_vaccine_parameters(
      safir_parameters = parameters,
      vaccine_ab_parameters = vax_pars,
      vaccine_set = vaccine_set,
      dose_period = dose_period,
      strategy_matrix = vaccine_coverage_strategy,
      next_dose_priority_matrix = next_dose_priority
    )
    
    parameters$mu_ab_infection <- mu_ab_infection
    
    # make VFR reduction vector and attach
    vfr_time1 <- as.Date(x = vfr_time1, format = "%m/%d/%Y")
    vfr_time2 <- as.Date(x = vfr_time2, format = "%m/%d/%Y")
    
    stopifnot(vfr_time1 < tmax_date)
    stopifnot(vfr_time2 < tmax_date)
    
    vfr_time1_day <- as.integer(difftime(vfr_time1, R0_t0 - 1))
    vfr_time2_day <- as.integer(difftime(vfr_time2, R0_t0 - 1))
    
    vfr_vector <- variant_fold_reduction_vector(parameters = parameters, dt = dt, vfr = vfr, vfr_time_1 = vfr_time1_day, vfr_time_2 = vfr_time2_day)
    parameters <- make_immune_parameters(parameters = parameters, vfr = vfr_vector, mu_ab_infection = mu_ab_infection, std10_infection = std10_infection)
    
    
    # apply modification to rate of hospitalisation and ICU
    
    v <- parameters$prob_hosp
    parameters$prob_hosp <- matrix(v, nrow=time_period, ncol=length(v), byrow=TRUE)
    mult <- c(rep(1,vfr_time1_day),seq(1,hosp_scal_omicron,length.out=vfr_time2_day-vfr_time1_day),rep(hosp_scal_omicron,time_period-vfr_time2_day))
    new <- parameters$prob_hosp*mult
    parameters$prob_hosp <- t(new)
    
    v <- parameters$prob_severe
    parameters$prob_severe <- matrix(v, nrow=time_period, ncol=length(v), byrow=TRUE)
    mult <- c(rep(1,vfr_time1_day),seq(1,ICU_scal_omicron,length.out=vfr_time2_day-vfr_time1_day),rep(ICU_scal_omicron,time_period-vfr_time2_day))
    new <- parameters$prob_severe*mult
    parameters$prob_severe <- t(new)
    
    
    # read in decay rate vector
    
    dr_vec <- readRDS(paste0("data/dr_vec_", fit, ".rds"))
    dr_vec_vaccine <- dr_vec %>%
      select(-t,-N)
    
    # assume that the decay rate for natural infection is the same as for D3 
    dr_vec_inf <- dr_vec %>%
      select(-t,-D1,-D2,-D3) 
    
    dr_vec_doses_m <- data.matrix(dr_vec_vaccine)
    dr_vec_inf_m <- data.matrix(dr_vec_inf)
    parameters <- make_independent_vaccine_infection_nat_parameters(parameters=parameters, dr_vec_doses=dr_vec_doses_m, dr_vec_inf=dr_vec_inf_m,max_ab_inf=5)
    
    ######################################################
    # run the simulation
    ######################################################
    
    # create variables
    timesteps <- parameters$time_period/dt
    
    # creates the categorical states and ages for the simulated population
    variables <- create_variables(pop = pop, parameters = parameters)
    variables <- create_vaccine_variables(variables = variables, parameters = parameters)
    variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)
    variables <- create_independent_nat_variables(variables = variables, parameters = parameters)
    
    # creates the list of events and attaches listeners which handle state changes and queue future events
    events <- create_events(parameters = parameters)
    events <- create_events_vaccination(events = events, parameters = parameters)
    attach_event_listeners(variables = variables, events = events, parameters = parameters, dt = dt)
    attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)
    attach_event_listeners_independent_nat(variables = variables, events = events, parameters = parameters, dt = dt)
    
    # renderer object is made
    renderer <- individual::Render$new(parameters$time_period)
    vaxx_renderer <- individual::Render$new(parameters$time_period)
    inf_renderer <- individual::Render$new(parameters$time_period)
    incidence_renderer <- individual::Render$new(timesteps)
    attach_tracking_listener_incidence(events=events, renderer = incidence_renderer)
    nat_renderer <- individual::Render$new(parameters$time_period)
    nat_inf_renderer <- individual::Render$new(parameters$time_period)
    sp_renderer <- individual::Render$new(parameters$time_period)
    hosp_render <- create_hosp_renderers(parameters = parameters)
    
    # track incidence
    attach_hosp_listeners(renderers = hosp_render, events = events)
    
    double_count_render_process_daily <- function(renderer, variable, dt) {
      stopifnot(inherits(variable, "DoubleVariable"))
      stopifnot(inherits(renderer, "Render"))
      function(t) {
        if ((t * dt) %% 1 == 0) {
          day <- as.integer(t * dt)
          nat <- exp(variable$get_values())
          quantiles_nat <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "lower", value = quantiles_nat[[1]], timestep = day)
          renderer$render(name = "median", value = quantiles_nat[[2]], timestep = day)
          renderer$render(name = "upper", value = quantiles_nat[[3]], timestep = day)
          renderer$render(name = "mean", value = mean(x = nat), timestep = day)
          }
      }
    }
    
    sp_render_process_daily <- function(renderer, variable1, variable2, dt) {
      stopifnot(inherits(variable1, "DoubleVariable"))
      stopifnot(inherits(variable2, "DoubleVariable"))
      stopifnot(inherits(renderer, "Render"))
      function(t) {
        if ((t * dt) %% 1 == 0) {
          day <- as.integer(t * dt)
          nat <- exp(variable1$get_values()) + exp(variable2$get_values())
          quantiles <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "nat_lower", value = quantiles[[1]], timestep = day)
          renderer$render(name = "nat_median", value = quantiles[[2]], timestep = day)
          renderer$render(name = "nat_upper", value = quantiles[[3]], timestep = day)
          renderer$render(name = "nat_mean", value = mean(x = nat), timestep = day)
          sp <-ifelse(nat>=ab_50,1,0)
          quantiles <- quantile(x = sp, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "sp_lower", value = quantiles[[1]], timestep = day)
          renderer$render(name = "sp_median", value = quantiles[[2]], timestep = day)
          renderer$render(name = "sp_upper", value = quantiles[[3]], timestep = day)
          renderer$render(name = "sp_mean", value = mean(x = sp), timestep = day)
          
        }
      }
    }
    
    # processes
    processes <- list(
      #    natural_immunity_ab_titre_process(parameters = parameters,variables = variables, dt = dt),
      independent_ab_titre_process(parameters=parameters, variables = variables, dt = dt),
      vaccination_process(parameters = parameters,variables = variables,events = events, dt = dt),
      infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events, dt = dt),
      categorical_count_renderer_process_daily(renderer = renderer, variable = variables$states, categories = variables$states$get_categories(),dt = dt),
      double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre, dt = dt ),
      double_count_render_process_daily(renderer = nat_inf_renderer, variable = variables$ab_titre_inf, dt = dt ),
      sp_render_process_daily(renderer = sp_renderer, variable1 = variables$ab_titre, variable2 = variables$ab_titre_inf, dt = dt),
      integer_count_render_process_daily(renderer = vaxx_renderer, variable = variables$dose_num, margin = 1:4, dt = dt))
    
    setup_events(parameters = parameters, events=events, variables = variables, dt=dt)
    
    simulation_loop_safir(
      variables = variables,
      events = events,
      processes = processes,
      timesteps = timesteps,
      variables_dont_update = c("discrete_age", "phase"),
      progress = TRUE
    )
    
     
    df <- renderer$to_dataframe()
    df_vacc <- vaxx_renderer$to_dataframe()
    df_inc <- incidence_renderer$to_dataframe()
    df_nat <- nat_renderer$to_dataframe() %>%
      rename(vax_ab_mean = mean,
             vax_ab_median = median,
             vax_ab_lower = lower,
             vax_ab_upper = upper)
    df_nat_inf <- nat_inf_renderer$to_dataframe() %>%
      rename(nat_ab_mean = mean,
             nat_ab_median = median,
             nat_ab_lower = lower,
             nat_ab_upper = upper)
    df_sp <- sp_renderer$to_dataframe()
    
    df_inc <- df_inc %>%
      mutate(timestep = floor((timestep-1)*dt)) %>% 
      mutate(incidence = coalesce(incidence,0))%>%
      group_by(timestep) %>%
      summarise(incidence = sum(incidence))
    
    df_hosp <- process_hosp_renderers(renderers = hosp_render, parameters = parameters)
    
    df_hosp <- df_hosp %>%
              mutate(hosp = hosp_get_live + hosp_get_die,
                     hosp_all = hosp + hosp_not_get_live + hosp_not_get_die,
                     ICU = ICU_get_live + ICU_get_die,
                     ICU_all = ICU + ICU_not_get_live + ICU_not_get_die,
                     timestep = day-1) %>%
              select(timestep,hosp,hosp_all,ICU,ICU_all)
    
    df <- left_join(df, df_vacc, by = c("timestep"))
    df <- left_join(df, df_inc, by = c("timestep"))
    df <- left_join(df, df_nat, by = c("timestep"))
    df <- left_join(df, df_nat_inf, by = c("timestep"))
    df <- left_join(df, df_sp, by = c("timestep"))
    df_rt <- as.data.frame(rt_out) %>%
      rename("timestep" = "Rt_tt")
    df <- left_join(df, df_rt, by = c("timestep"))
    df <- left_join(df, df_hosp, by = c("timestep"))
    
    
    
    # summarise
    saf_reps_summarise <- df %>%
      mutate(IMild_count = IMild_count + IAsymp_count) %>%
      select(-IAsymp_count) %>%
      pivot_longer(cols = contains(c("count", "Rt","incidence","nat","sp", "vax","hosp","ICU")), names_to = "compartment") %>%
      
      filter(compartment %in% c("D_count", "X1_count", "X2_count", "X3_count", "R_count", "IMild_count", "ICase_count", "Rt", "incidence", 
                                "vax_ab_mean", "vax_ab_lower", "vax_ab_upper", "nat_ab_mean", "nat_ab_lower", "nat_ab_upper",
                                "nat_mean", "nat_lower","nat_upper","sp_mean","sp_lower", "sp_upper", 
                                "hosp", "hosp_all", "ICU", "ICU_all" )) %>%
      group_by(compartment) %>%
      mutate(value = if_else(compartment == "D_count", value - lag(value), value),
             value = if_else(is.na(value), 0, value)) %>%
      ungroup() %>%
      pivot_wider(id_cols = timestep, names_from = "compartment", values_from = "value")  %>%
      mutate(deaths = sum(D_count[timestep >= days_to_vacc_start]),
             cum_hosp = sum(hosp[timestep >= days_to_vacc_start]),
             cum_hosp_all = sum(hosp_all[timestep >= days_to_vacc_start]),
             cum_ICU = sum(ICU[timestep >= days_to_vacc_start]),
             cum_ICU_all = sum(ICU_all[timestep >= days_to_vacc_start]),
             inc = sum(incidence[timestep >= days_to_vacc_start]),
             doses = X1_count + X2_count * 2 + X3_count * 3,
             total_doses = max(doses)) %>%
      ungroup() %>%
      nest(cols = c(timestep, D_count, X1_count, X2_count, X3_count, R_count, IMild_count, ICase_count, doses, Rt, incidence,
                    vax_ab_mean, vax_ab_lower, vax_ab_upper, nat_ab_mean, nat_ab_lower, nat_ab_upper,
                    nat_mean, nat_lower, nat_upper, sp_mean, sp_lower, sp_upper,
                    hosp, hosp_all, ICU, ICU_all)) %>%
      mutate(scenario = scenario)
    
    # get prop seropositive when vaccination starts
    prop_R <- df %>%
      select(timestep, sp_mean) %>%
      filter(timestep == days_to_vacc_start) %>%
      mutate(prop_R = round(sp_mean * 100,2)) %>%
      select(prop_R) %>%
      mutate(scenario = scenario)
    
    saf_reps_summarise <- left_join(saf_reps_summarise, prop_R, by = "scenario")
    
    # Save output
    # to create folder:
    output_address <- dir.create(paste0("raw_outputs/output_", name), recursive = TRUE)
    saveRDS(saf_reps_summarise, paste0(paste0("raw_outputs/output_", name, "/scenario_", scenario),".rds"))
    
    #if folder present:
    #output_address <- paste0("raw_outputs/output_", name, "/scenario_", scenario, ".rds")
    #saveRDS(saf_reps_summarise, output_address)
  }


name <- "india_run_1"

###############################################################################################################
## we directly read in vaccine efficacy and decay parameters from fits made to the English VE data
## fit1 gives the path for the posterior parameter set
## fit is the label for the decay rate in NATs that has been modified to be slower following dose 3
## and for natural immunity. it is stored in an RDS dataset. 
###############################################################################################################

fit1 <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"
fit <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_SD1"

###############################################################################################################
## the code is now adapted to read in an Rt profile - so that we can link to nimue fits
## note that the safir interpolation is currently active and throws an error for daily estimates (issue flagged)
## but we don't need interpolation if we have daily estimates - just need to extend to have Rt from the start time
## through to the end time of the simulation
###############################################################################################################

R_profile <- read_csv("data/category_2_Rt.csv") 
saveRDS(R_profile,"data/R_profile_india.rds")


#### Get vaccine parameters  #############################################################
#### note: will need to change the booster dose as using AZ not PF #######################

vaccine <- "Oxford-AstraZeneca"
vacc_names <- data.frame(vaccine = c("Pfizer", "Oxford-AstraZeneca"), vacc = c("PF", "AZ"))
vaccine_n <- vaccine
vacc_params <- readRDS(paste0("data/param_list_",fit1,".rds")) %>%
  rename(vacc = vaccine) %>%
  left_join(vacc_names, by = "vacc") %>%
  filter(vaccine == vaccine_n) %>%
  mutate(std10 = 0.44) %>%
  filter(vfr > 1) %>%
  select(-c(vacc))

#define these for use when testing code run function
test <- 0
if(test){
  source("R/utils.R")
  source("R/vaccine_strategy.R")
  source("R/run_function_india.R")
  mu_ab_d1 <- vacc_params$mu_ab_d1[1]
  mu_ab_d2 <- vacc_params$mu_ab_d2[1]
  k <- vacc_params$k[1]
  dose_3_fold_increase <- vacc_params$dose_3_fold_increase[1]
  hl_s <- vacc_params$hl_s[1]
  hl_l <- vacc_params$hl_l[1]
  period_s <- vacc_params$period_s[1]
  period_l <- 365
  ab_50 <- vacc_params$ab_50[1]
  ab_50_severe <- vacc_params$ab_50_severe[1]
  std10 <- vacc_params$std10[1]
  std10_infection <- std10
}

#### Set up other simulation parameters  ##############################################

target_pop <- 5e5
dt <- 0.25
repetition <- 1:5
seeding_cases <- 10
ab_model_infection <- TRUE
max_Rt <- 4
mu_ab_infection <- 1.7
hs_constraints <- "Present"

#### Omicron specific parameters ###############################################

max_Rt_omicron <- 4.25 #c(4.25,4.5,4.75)
vfr_time1 <- "11/27/2021"
vfr_time2 <- "12/31/2021"
vfr <- 3.9 #unique(vacc_params$vfr)
central_vfr <- vfr[1]
hosp_scal_omicron <- 0.5
ICU_scal_omicron <-  0.42

### vaccine roll-out parameters ################################################

vacc_start <- "4/1/2021"
vaccine_doses <- 2 #c(2,3)
max_coverage <- 0.8
age_groups_covered <- 9 #c(5, 9)
vacc_per_week <- 0.02 #c(0.02, 0.015)
t_d3 <- 180

#### Create scenarios ##########################################################

scenarios <- expand_grid(fit = fit,
                         hs_constraints = hs_constraints,
                         target_pop = target_pop,
                         vaccine_doses = vaccine_doses,
                         vaccine = vaccine,
                         max_coverage = max_coverage,
                         age_groups_covered = age_groups_covered,
                         vacc_start = vacc_start,
                         dt = dt,
                         repetition = repetition,
                         seeding_cases = seeding_cases,
                         vacc_per_week = vacc_per_week,
                         ab_model_infection = ab_model_infection,
                         t_d3 = t_d3,
                         max_Rt = max_Rt,
                         max_Rt_omicron = max_Rt_omicron,
                         vfr = vfr,
                         vfr_time1 = vfr_time1,
                         vfr_time2 = vfr_time2,
                         mu_ab_infection = mu_ab_infection,
                         hosp_scal_omicron = hosp_scal_omicron,
                         ICU_scal_omicron = ICU_scal_omicron) %>%
unique()

scenarios$scenario <- 1:nrow(scenarios)
scenarios$name <- name
# add vaccine parameter list
scenarios <- left_join(scenarios, vacc_params, by = c("vaccine", "vfr")) 

nrow(scenarios)

write_csv(scenarios, paste0("scenarios_", name, ".csv"))

## test on PC

source("R/run_function_india.R")
source("R/utils.R")
source("R/vaccine_strategy.R")
plan(multicore, workers = 4)
system.time({out <- future_pmap(scenarios, run_scenario, .progress = TRUE)})

#### Run the model on cluster ###############################################
# Load functions
sources <- c("R/run_function_abmodel_omicron_india.R", "R/utils.R", "R/vaccine_strategy.R")
src <- conan::conan_sources(c("mrc-ide/safir", "mrc-ide/squire", "mrc-ide/nimue"))
ctx <- context::context_save("context",
                             sources = sources,
                             packages = c("tibble", "dplyr", "tidyr", "countrycode", "safir", "nimue", "squire", "data.table"),
                             package_sources = src)

config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--didemrchnb")
#config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--dideclusthn")

# Create the queue
run <- didehpc::queue_didehpc(ctx, config = config)
# Summary of all available clusters
# run$cluster_load(nodes = FALSE)
# Run
runs <- run$enqueue_bulk(scenarios, run_scenario, do_call = TRUE, progress = TRUE)
runs$status()

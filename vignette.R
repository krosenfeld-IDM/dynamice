# library (dynamice)
devtools::load_all()
library (ggplot2)
library (scales)
library (foreach)
library (iterators)
library (parallel)
library (doParallel)
library (countrycode)
library (ggpubr)
library (lhs)
library (truncnorm)

# include functions (delete when package is ready)
source ("R/logs.R")
source ("R/utils.R")
source ("R/functions.R")
source ("R/dxplot.R")

# load data (delete when package is ready)
load(file = "data/data_pop.rda")
load(file = "data/data_cfr_portnoy.rda")
load(file = "data/data_cfr_wolfson.rda")
load(file = "data/data_contact_polymod.rda")
load(file = "data/data_contact_syn.rda")
load(file = "data/data_r0.rda")
load(file = "data/data_timeliness.rda")
load(file = "data/data_lexp_remain.rda")
load(file = "data/data_template.rda")

# assume perfect timeliness: all receive MCV1 at 9-month
adj.timeliness = TRUE  # FALSE
if (adj.timeliness){
  data_timeliness [!is.na(age) & age < 52*(9/12), timeliness := 0]
  data_timeliness [!is.na(age) & age >= 52*(9/12), timeliness := 1]
}
# assume a fixed R0
adj.fixR0 = 16  # NA
if (!is.na(adj.fixR0)){
  data_r0 [ , r0 := adj.fixR0]
}

# setup the scenarios

fortran_names <- c("vaccine2019_sia_adjBeta_mcv2_rSIAdpd.exe",          # most updated version
                   "vaccine2019_sia_adjBeta_mcv2_nodpdSIA.exe",
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAdpd.exe",
                   "vaccine2019_sia_adjBeta_mcv2_nodpdSIA_stepVE.exe",
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAdpdSA.exe", # sensitivity analysis for SIA delivery
                   "vaccine2019_sia_adjBeta_mcv2_stepVE_rSIAnever.exe"  # alternative assumption for SIA delivery
)

contact_mats <- c("syn", "polymod", "prpmix", "unimix")

main <- list (
  scenarios    = c("Base",          # (1) Base: step change VE + random SIA delivery to zero-dose children
                   "Efficacy",      # (2) Base + linear VE
                   "SIAdeliver",    # (3) Base + non-random SIA delivery to zero-dose children
                   "ContactSyn",    # (4) Base + synthetic contact patterns
                   "ContactPrp",    # (5) Base + proportional mixing
                   "ContactUni",    # (6) Base + uniform mixing
                   "Update",        # (7) Update: linear VE + synthetic matrix + non-random SIA delivery
                   "SIAdeliverSA",  # (8) Base + sensitivity analysis for non-random SIA delivery
                   "SIAnever"       # (9) Base + SIA delivery assuming 7.7% never reached
  ),
  fortran_model = fortran_names [c(4, 2, 3, 4, 4, 4, 1, 5, 6)],
  contact_mat   = contact_mats  [c(2, 2, 2, 1, 3, 4, 1, 2, 2)],
  step_ve       =                c(T, F, T, T, T, T, F, T, T)
)

# evaluate the scenario
imain <- 7
print (paste0 ("Evaluation started: ", main$scenarios[imain]))

var <- list (
  # vaccine coverage
  vaccine_coverage_folder           = "vac_coverage_top10/base case/",
  coverage_prefix                   = "coverage",
  touchstone                        = "_",
  antigen                           = NULL,
  vaccine_coverage_subfolder        = "scenarios/",

  # disease burden
  central_burden_estimate_folder    = "central_burden_estimate/",
  stochastic_burden_estimate_folder = "stochastic_burden_estimate/",

  # diagnostic plots folder
  plot_folder                       = "plots/",

  # modelling group name
  group_name                        = NULL,

  # log file name
  log_name                          = "test_log_1",

  # countries - specify iso3 codes to analyse only these countries
  #             or set it to "all" to analyse all included countries
  countries                         = c("PAK"),#

  cluster_cores                     = 1,    # number of cores
  psa                               = 0     # psa runs; 0 for single central run
)

# for central run: set number of runs to 0 (var$psa = 0)
# for stochastic runs: set number of runs to greater than 0 (eg: var$psa = 200)
# burden_estimate_folder (central or stochastic)
if (var$psa == 0) {
  burden_estimate_folder <- var$central_burden_estimate_folder
} else {
  burden_estimate_folder <- var$stochastic_burden_estimate_folder
}

# folders for burden estimate results
dir.create (file.path (paste0 (getwd(), "/", burden_estimate_folder, "Wolfson/")), recursive = T)

dir.create (file.path (paste0 (getwd(), "/", burden_estimate_folder, "Portnoy/")), recursive = T)


vaccine_strategies <- c("no-vacc",    # (1) no vaccination (set vaccination and using_sia to 0)
                        "mcv1",       # (2) MCV1 only
                        "mcv2",       # (3) MCV1 & MCV2
                        "sia-conti",  # (4) MCV1 & MCV2 and SIAs
                        "sia-stop"    # (5) MCV1 & MCV2 and SIAs (stop SIAs after 2019)
)

# corresponding vaccination strategies used for Portnoy's CFR
portnoy_scenarios <- c("no-vaccination",    # (1) no vaccination (set vaccination and using_sia to 0)
                       "mcv1-default",      # (2) MCV1 only
                       "mcv2-default",      # (3) MCV1 & MCV2
                       "campaign-default",  # (4) MCV1 & MCV2 and SIAs
                       "campaign-default"   # (5) MCV1 & MCV2 and SIAs (stop SIAs after 2019)
)

# set SIAs and vaccination parameters for each scenario to minimize errors for running
set_sia         <- c (0, 0, 0, 1, 1)
set_vaccination <- c (0, 1, 2, 2, 2)

# Prepare coverage
update_coverage_inputs = FALSE
if (update_coverage_inputs) {
  # (1) create vaccine coverage file (0% coverage) for no vaccination scenario
  create_no_vaccination_coverage_file (
    no_vaccination_coverage_file = paste0 (var$vaccine_coverage_folder,
                                           "coverage_", scenarios[1],".csv"),
    vaccination_coverage_file    = paste0 (var$vaccine_coverage_folder,
                                           "coverage_", scenarios[2],".csv")
  )

  # (2) create a vaccine coverage file for no future SIA scenario
  continuous_sia_coverage <- fread (paste0(var$vaccine_coverage_folder,
                                           "coverage_", scenarios[4],".csv"))
  stop_sia_coverage <- continuous_sia_coverage [!(year > 2019 & activity_type == "campaign")]
  fwrite (stop_sia_coverage, file = paste0(var$vaccine_coverage_folder,
                                           "coverage_", scenarios[5],".csv"))

  # (3) generate 2 vaccine coverage files per scenario for routine and SIA
  for (index in 1:5) {
    create_vaccine_coverage_routine_sia (
      vaccine_coverage_folder    = var$vaccine_coverage_folder,
      vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
      coverage_prefix            = var$coverage_prefix,
      touchstone                 = var$touchstone,
      antigen                    = var$antigen,
      scenario_name              = scenarios [index],
      rev_cov                    = FALSE # coverage data have been pre-processed
    )
  }
}

# create samples for probabilistic analysis
if (var$psa > 0) {
  psadat <- create_PSA_data (psa             = var$psa,
                             seed_state      = 1,
                             psadat_filename = "psa_variables.csv")
}

# ger results
# selected scenario
index <- 4
scenario_name  <- vaccine_strategies [index]
print (scenario_name)

scenario_number <- sprintf ("scenario%02d", index)

# copy data files
data_timeliness = copy(data_timeliness)
data_r0 = copy(data_r0)
data_pop = copy(data_pop)
data_template = copy(data_template)


# run model and estimate cases
burden_estimate_file <- runScenario (
  vaccine_coverage_folder    = var$vaccine_coverage_folder,
  vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
  coverage_prefix            = var$coverage_prefix,
  antigen                    = var$antigen,
  scenario_name              = scenario_name,
  save_scenario              = scenario_number,
  burden_estimate_folder     = burden_estimate_folder,
  group_name                 = var$group_name,
  log_name                   = var$log_name,
  countries                  = var$countries,
  cluster_cores              = var$cluster_cores,
  psa                        = var$psa,
  vaccination                = set_vaccination [index],
  using_sia                  = set_sia         [index],
  measles_model              = main$fortran_model [imain],
  debug_model                = TRUE,
  contact_mat                = main$contact_mat [imain],
  step_ve                    = main$step_ve     [imain],
  sim_years                  = 1980:2100
)

# estimate deaths and DALYs by Wolfson and Portnoy's methods
estimate_deaths_dalys (cfr_option             = "Wolfson",
                       burden_estimate_file   = burden_estimate_file,
                       burden_estimate_folder = burden_estimate_folder,
                       psa                    = var$psa)

estimate_deaths_dalys (cfr_option             = "Portnoy",
                       burden_estimate_file   = burden_estimate_file,
                       burden_estimate_folder = burden_estimate_folder,
                       vimc_scenario          = portnoy_scenarios[index],
                       portnoy_scenario       = "s6",  # constant CFR after 2018
                       psa                    = var$psa)

# Get diagnostic plots
if (var$psa == 0) {

  dir.create (file.path (paste0 (getwd(), "/", var$plot_folder)), recursive = T)

  diagnostic_plots (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    coverage_prefix            = var$coverage_prefix,
    touchstone                 = var$touchstone,
    antigen                    = var$antigen,
    scenarios                  = vaccine_strategies[index],
    base_scenario              = "no-vaccination",
    burden_estimate_folder     = var$central_burden_estimate_folder,
    plot_folder                = var$plot_folder,
    group_name                 = var$group_name,
    countries                  = var$countries,
    cfr_options                = c("Wolfson", "Portnoy"),
    psa                        = var$psa,
    start_year                 = 2000,
    end_year                   = 2050,
    compare_plots              = FALSE
  )
}

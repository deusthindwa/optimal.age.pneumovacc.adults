# load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", 
                        "here",
                        "scales", 
                        "magrittr",
                        "mvtnorm"))

options(stringsAsFactors = FALSE)
setwd(here::here())

source(here::here("script", "load_data.R"))

# model incidences in each age group using 
source(here::here("script", "incidence.R"))

# estimate vaccine impact against all IPD serotypes
source(here::here("script", "metacurve.R"))

# generate scenarios by serogroup, country, VE, age and waning
source(here::here("script", "vaccination_scenarios.R"))

# compute vaccine impact
source(here::here("script", "vaccine_impact.R"))

source(here::here("script", "plots.R"))
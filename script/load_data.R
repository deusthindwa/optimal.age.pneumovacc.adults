# load the IPD cases
ipd <- readr::read_csv(here("data", "EW_ipd_incid.csv"))
source(here::here("script", "pops.R"))
ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup,1,2)))

Nsims <- 1e3 # number of simulations to use for all uncertainty analysis
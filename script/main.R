# load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", 
                        "here",
                        "scales", 
                        "magrittr"))
setwd(here::here())


# load the IPD cases
ipd <- readr::read_csv(here("data", "EW_ipd_incid.csv"))
source(here::here("script", "pops.R"))
ipd <- dplyr::mutate(ipd, agey = readr::parse_number(substr(agegroup,1,2)))


# estimate the rest of parameters using a simple linear model
theta0 <- min(ipd$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=ipd)  
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]


# fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
# we parameterise in terms of log-rates to ensure the estimates are positive
ipd_models <- ipd %>% 
  split(.$serogroup) %>%
  purrr::map(~nls(data = .x, 
                  incidence ~ exp(log_alpha) * exp(exp(log_beta)*agey) + exp(log_theta),
                  nls.control(maxiter = 200),
                  start = list(log_alpha = log(alpha0),
                               log_beta  = log(beta0),
                               log_theta = log(theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))

ipd_curves <- ipd_models %>%
  purrr::map_df(~dplyr::mutate(ipd_x,
                               incidence = predict(object = .x,
                                                   newdata = ipd_x)),
                .id = "serogroup")


# calculate scaled incidence
ipd_scaled <- ipd %>% 
  dplyr::group_by(serogroup) %>%
  dplyr::mutate(p = incidence/sum(incidence))

# generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- dplyr::inner_join(ipd_curves, countries_df, by = "agey") %>%
  dplyr::filter(serogroup != "All serotypes") %>%
  dplyr::mutate(cases = incidence/1e5*ntotal)

# estimate vaccine impact against all IPD serotypes
source(here::here("script", "metacurve.R"))

VE_table <- tibble::add_row(VE_table, 
                            Study       = "None", 
                            rate        = 0, 
                            `Half-life` = Inf)

# make a list check it twice
initial_VE <- function(age, serogroup, age_dep = FALSE, scale = 1){
  # age-dependent vaccine efficacy at time of vaccination
  # may be superseded
  dplyr::case_when(serogroup == "PCV13" ~ 
                     dplyr::case_when(age_dep == FALSE ~ NA_real_,
                                      age >= 75 ~ 0.30,
                                      age >= 65 ~ 0.36,
                                      age >= 55 ~ 0.54,
                                      TRUE    ~ NA_real_),
                   # if no age dependency, we need to just use the value from
                   # the relevant study, which we can handle outside this
                   serogroup == "PPV23" ~
                     dplyr::case_when(age_dep == FALSE ~ NA_real_,
                                      age >= 75 ~ 0.30,
                                      age >= 65 ~ 0.36,
                                      age >= 55 ~ 0.54,
                                      TRUE    ~ NA_real_),
                   TRUE                 ~ NA_real_)/scale
}

# generate scenarios by serogroup, country, VE, age and waning
options(stringsAsFactors = FALSE)

scenarios <- list(`1` = data.frame(Study.waning = "Andrews (2012)",
                                   Study.VE     = "Andrews (2012)"),
                  `2` = data.frame(Study.waning = "Djennad (2018)",
                                   Study.VE     = "Djennad (2018)"),
                  `3` = data.frame(Study.waning = "Andrews (2012)"),
                  `4` = data.frame(Study.waning = "Djennad (2018)")) %>%
  dplyr::bind_rows(.id = "scenario") %>%
  dplyr::mutate(age_dep = is.na(Study.VE)) %>%
  tidyr::crossing(Vac.age = seq(55, 85, by = 5),
                  serogroup = c("PPV23",
                                "PCV13")) %>%
  dplyr::mutate(Study.waning = ifelse(serogroup == "PCV13" & scenario %in% c(2,4),
                                      "Andrews (2012)",
                                      Study.waning)) %>%
  dplyr::left_join(
    dplyr::select(VE_table,
                  Study.waning = Study,
                  rate),
    by = "Study.waning") %>%
  dplyr::left_join(
    dplyr::select(VE_table,
                  Study.VE = Study,
                  VE),
    by = "Study.VE") %>%
  dplyr::mutate(age_dep = is.na(Study.VE)) %>% 
  dplyr::mutate(
    VE    = ifelse(test = age_dep,
                   yes  = initial_VE(Vac.age, serogroup, age_dep),
                   no   = VE),
    rate  = ifelse(test = serogroup == "PCV13" & scenario %in% c(1,3),
                   yes  = 0,
                   no   = rate),
    delay = ifelse(test = serogroup == "PCV13" & scenario %in% c(2,4),
                   yes  = 5,
                   no   = 0))

# we want the curve to be at VE if age >= vac.age + delay
# when delay > 0, we want to subtract delay off
VE_by_Vac.age <- 
  scenarios %>%
  dplyr::inner_join(Cases) %>%
  dplyr::mutate(agey_since = agey - Vac.age) %>%
  dplyr::mutate(agey_since = ifelse(test = delay == 0,
                                    yes  = agey_since,
                                    no   = pmax(0, agey_since - delay)
  )) %>% 
  dplyr::mutate(Vaccine_Efficacy = VE*exp(rate*(1 + agey_since))) %>%
  dplyr::mutate(value = ifelse(agey < Vac.age, 0, Vaccine_Efficacy)) %>%
  dplyr::mutate(Impact = value*cases) 


# compute vaccine impact
impact_by_age_to_plot <- 
  VE_by_Vac.age %>%
  dplyr::group_by(serogroup,
                  Study.waning,
                  Study.VE,
                  age_dep, delay,
                  rate,
                  Vac.age, Country,
                  scenario) %>%
  dplyr::summarise(Impact = sum(Impact)) %>%
  dplyr::mutate(Waning = dplyr::case_when(
    rate         == 0 ~ "No waning",
    Study.waning == "Andrews (2012)" ~ "Fast waning",
    Study.waning == "Djennad (2018)" ~ "Slow waning"),
    Waning = ifelse(delay > 0,
                    paste(Waning, sprintf("\n(%i years' delay)", delay)),
                    Waning)) 


# compute maximum vaccine impact
impact_by_age_to_plot_max <- 
  impact_by_age_to_plot %>%
  dplyr::group_by_at(.vars = dplyr::vars(-c(Vac.age, Impact))) %>%
  dplyr::filter(Impact == max(Impact))


# vaccine impact per 10000 older adults
Area <- rename(subset(countries_df, select = c("Country","agey","ntotal")), c("Vac.age"="agey"))
impact_validated <- merge(impact_by_age_to_plot, Area, by=c("Country","Vac.age"))


# 65y old programme impact (%) at 70% coverage
A65 <- subset(impact_by_age_to_plot, serogroup=="PPV23" & Country == "England/Wales")
(A65[A65$Waning=="Fast waning" & A65$scenario==3 & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Fast waning" & A65$scenario==3,]$Impact))*.7
(A65[A65$Waning=="Fast waning" & A65$scenario==1 & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Fast waning" & A65$scenario==1,]$Impact))*.7
(A65[A65$Waning=="Slow waning" & A65$scenario==4 & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Slow waning" & A65$scenario==4,]$Impact))*.7
(A65[A65$Waning=="Slow waning" & A65$scenario==2 & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Slow waning" & A65$scenario==2,]$Impact))*.7



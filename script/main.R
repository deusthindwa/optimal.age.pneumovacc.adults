
# load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", 
                        "here",
                        "scales", 
                        "magrittr",
                        "mvtnorm",
                        "zoo",
                        "patchwork"))

<<<<<<< HEAD
#bootstrap to obtain confidence intervals for fitted values

uncertainty <- function(z){
  agey <- ipd_curves$agey[z]
  incidence <- ipd_curves$incidence[z]
  ipd_modelx <- nls(data =ipd_curves, 
            incidence ~ exp(log_alpha) * exp(exp(log_beta)*agey) + exp(log_theta), 
            start = list(log_alpha = log(0.0409), log_beta = log(0.0772), log_theta = log(1.71)),
            nls.control(maxiter = 200))
  
  predict(ipd_modelx, newdata = list(agey = agey))
}

intervals <-as_data_frame(replicate(1000, uncertainty(sample.int(108, replace = TRUE))))
interval <- apply(intervals, 1, quantile, probs = c(0.05, 0.95))
l_int <- apply(intervals, 1, quantile, probs = 0.05)
u_int <- apply(intervals, 1, quantile, probs = 0.95)


summary(intervals)











=======
options(stringsAsFactors = FALSE)
setwd(here::here())
>>>>>>> c571dba80c4d93d1981f10e7603ffc0586f595bd

source(here::here("script", "load_data.R"))

# model incidences in each age group using 
source(here::here("script", "incidence.R"))

# estimate vaccine impact against all IPD serotypes
source(here::here("script", "metacurve.R"))

# generate scenarios by serogroup, country, VE, age and waning
source(here::here("script", "vaccination_scenarios.R"))

# compute vaccine impact
<<<<<<< HEAD
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
(A65[A65$Waning=="Fast waning" & A65$age_dep==TRUE & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Fast waning" & A65$Vac.age==65,]$Impact))*.7
(A65[A65$Waning=="Fast waning" & A65$age_dep==FALSE & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Fast waning" & A65$Vac.age==65,]$Impact))*.7
(A65[A65$Waning=="Slow waning" & A65$age_dep==TRUE & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Slow waning" & A65$Vac.age==65,]$Impact))*.7
(A65[A65$Waning=="Slow waning" & A65$age_dep==FALSE & A65$Vac.age==65,]$Impact/sum(A65[A65$Waning=="Slow waning" & A65$Vac.age==65,]$Impact))*.7

=======
source(here::here("script", "vaccine_impact.R"))
>>>>>>> c571dba80c4d93d1981f10e7603ffc0586f595bd

# make plots
source(here::here("script", "plots.R"))

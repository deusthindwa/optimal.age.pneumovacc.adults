
VE_table <- tibble::add_row(VE_table, 
                            Study       = "None", 
                            rate        = 0, 
                            `Half-life` = Inf)


# make a list check it twice
initial_VE <- function(age, serogroup, age_dep = FALSE){
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
                     TRUE                 ~ NA_real_)
}

df_from_study_ <- distinct(df_from_study, Study, VE, rate, sim)


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
    crossing(sim = 1:Nsims) %>%
    dplyr::left_join(
        dplyr::select(df_from_study_,
                      Study.waning = Study,
                      sim,
                      rate),
        by = c("sim","Study.waning")) %>%
    dplyr::left_join(
        dplyr::select(df_from_study_,
                      Study.VE = Study,
                      VE,
                      sim),
        by = c("sim","Study.VE")) %>%
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

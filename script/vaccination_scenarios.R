
VE_table <- tibble::add_row(VE_table, 
                            Study       = "None", 
                            rate        = 0, 
                            `Half-life` = Inf)


# make a list check it twice
initial_VE <- function(age, serogroup, age_dep = FALSE){
    # age-dependent vaccine efficacy at time of vaccination
    # may be superseded
    dplyr::case_when(serogroup == "PCV13" ~ 
                                    dplyr::case_when(age_dep == FALSE ~ 1,
                                                     age >= 75 ~ 0.30/0.54,
                                                     age >= 65 ~ 0.36/0.54,
                                                     age >= 55 ~ 0.54/0.54,
                                                     TRUE    ~ NA_real_),
                                # if no age dependency, we need to just use the value from
                                # the relevant study, which we can handle outside this
                                serogroup == "PPV23" ~
                                    dplyr::case_when(age_dep == FALSE ~ 1,
                                                     age >= 75 ~ 0.30/0.54,
                                                     age >= 65 ~ 0.36/0.54,
                                                     age >= 55 ~ 0.54/0.54,
                                                     TRUE    ~ NA_real_),
                                TRUE                 ~ NA_real_)
}

df_from_study_ <- distinct(df_from_study, Study, VE, rate, sim)


# probably just need to do the single Study
# maybe just  double check that we're doing zero for 1,3 PCV
# this whole thing could just be an expand.grid
# 

scenarios <- list(`1` = data.frame(Study.waning = "Andrews (2012)",
                                   Study.VE     = "Andrews (2012)"),
                  `2` = data.frame(Study.waning = "Djennad (2018)",
                                   Study.VE     = "Djennad (2018)"),
                  `3` = data.frame(Study.waning = "Andrews (2012)",
                                   Study.VE     = NA),
                  `4` = data.frame(Study.waning = "Djennad (2018)",
                                   Study.VE     = NA)) %>%
    dplyr::bind_rows(.id = "scenario") %>%
    dplyr::mutate(age_dep = scenario >= 3) %>%
    tidyr::crossing(Vac.age = seq(55, 85, by = 5),
                    serogroup = c("PPV23",
                                  "PCV13")) 

scenarios <- expand.grid(
    data.frame(initial = c("Andrews (2012)", "Djennad (2018)"),
                                    serogroup = c("PCV13", "PPV23"),
                                    age_dep = c(TRUE, FALSE)),
                         stringsAsFactors = FALSE) %>%
    arrange(serogroup, age_dep) %>%
    mutate(waning = case_when(
        initial == "Andrews (2012)" ~ "Fast",
        serogroup == "PCV13" & initial == "Djennad (2018)" ~ "None",
        serogroup == "PPV23" & initial == "Djennad (2018)" ~ "Slow",
        TRUE ~ "Unknown")) %>%
    mutate(delay = ifelse(serogroup == "PCV13" & waning == "Fast",
                          5,
                          0),
           initial = ifelse(serogroup == "PCV13", "Andrews (2012)", initial)) %>%
    rename(Study.VE = initial) %>%
    mutate(Study.waning = case_when(waning == "Fast" ~ "Andrews (2012)",
                                    waning == "Slow" ~ "Djennad (2018)",
                                    waning == "None" ~ "None",
                                    TRUE ~ NA_character_)) %>%
    select(-waning)
    
scenarios %<>%
    crossing(sim = 1:Nsims,
             Vac.age = seq(55, 85, by = 5)) %>%
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
    dplyr::mutate(rate  = replace_na(data = rate, 0),
                  scale = initial_VE(Vac.age, serogroup, age_dep),
                  VE    = scale*VE) %>%
    select(-scale)
    

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


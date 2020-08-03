impact_by_age_to_plot <- 
    VE_by_Vac.age %>%
    dplyr::group_by(serogroup,
                    Study.waning,
                    Study.VE,
                    age_dep,
                    delay,
                    sim,
                    Vac.age,
                    Country) %>%
    dplyr::summarise(Impact = sum(Impact)) %>%
    dplyr::mutate(Waning = dplyr::case_when(
        Study.waning == "None"           ~ "No waning",
        Study.waning == "Andrews (2012)" ~ "Fast waning",
        Study.waning == "Djennad (2018)" ~ "Slow waning"),
        Waning = ifelse(delay > 0,
                        paste(Waning, sprintf("\n(%i years' delay)", delay)),
                        Waning)) 


# # compute maximum vaccine impact
# impact_by_age_to_plot_max <- 
#     impact_by_age_to_plot %>%
#     dplyr::group_by_at(.vars = dplyr::vars(-c(Vac.age, Impact))) %>%
#     dplyr::filter(Impact == max(Impact))


# vaccine impact per 10000 older adults
impact_validated <-
    dplyr::select(countries_df, Country, agey, ntotal) %>%
    dplyr::rename(Vac.age = agey) %>%
    dplyr::inner_join(impact_by_age_to_plot,
                      by=c("Country","Vac.age")) %>%
    mutate(Impact = Impact*10000/ntotal) %>%
    nest(data = c(sim, Impact)) %>%
    mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025,
                                                           0.5,
                                                           0.975)))) %>%
    unnest_wider(Q)



# 65y old programme impact (%) at 70% coverage
coverage <- 0.7

A65 <- impact_by_age_to_plot %>%
    dplyr::filter(serogroup == "PPV23",
                  Country   == "England/Wales",
                  Vac.age >= 65) %>% 
                  inner_join(subset(Cases, select = c("serogroup", "Vac.age", "Country", "sim", "cases")))

impact_65y_70pc <-
    dplyr::group_by(A65, Waning, sim) %>%
    dplyr::mutate(value = coverage*Impact/sum(cases),
                  Waning = sub(pattern     = "\\swaning", 
                               replacement = "", 
                               x           = Waning)#,
                  #Waning = sprintf("%s (%0.3f)", Waning, rate)
                  ) %>% 
    dplyr::filter(Vac.age == 65) %>%
    select(-delay, -Impact, -cases) %>%
    #dplyr::mutate(value = scales::percent(value, 0.1)) %>%
    nest(data = c(sim, value)) %>%
    mutate(Q = map(data, ~quantile(.x$value, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q) %>%
    ungroup %>%
    mutate_at(.vars = vars(contains("%")),
              .funs = ~scales::percent(., 0.1)) %>%
    select(Waning,
           `Age dep.` = age_dep,
           contains("%")) %>%
    group_by_at(.vars = vars(-contains("%"))) %>%
    transmute(Impact = sprintf("%s (%s, %s)", `50%`, `2.5%`, `97.5%`))


readr::write_csv(x    = impact_65y_70pc, 
                 path = here("output", "impact_65y_70pc.csv"))

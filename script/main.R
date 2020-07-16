#load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", "curl",
                        #"Hmisc",
                        "here",
                        "scales", "magrittr"))

setwd(here::here())

#load the IPD cases

ipd    <- read_csv(here("data", "EW_ipd_incid.csv"))

source(here("script", "pops.R"))

ipd <- mutate(ipd, agey = readr::parse_number(substr(agegroup,1,2)))

#estimate the rest of parameters using a simple linear model
theta0 <- min(ipd$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=ipd)  
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]

#fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
ipd_models <- ipd %>% 
  split(.$serogroup) %>%
  map(~nls(data = .x, 
           incidence ~ exp(log_alpha) * exp(exp(log_beta)*agey) + exp(log_theta),
           nls.control(maxiter = 200),
           start = list(log_alpha = log(alpha0),
                        log_beta  = log(beta0),
                        log_theta = log(theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))

ipd_curves <- ipd_models %>%
  map_df(~mutate(ipd_x, incidence = predict(object = .x, newdata = ipd_x)),
         .id = "serogroup")

#plot fitted curves with backward or forward extrapolation
incidence_plot <- ggplot(data = ipd,
                         aes(x = agey,
                             y = incidence,
                             color = serogroup)) + 
  geom_point(size=2.5) + 
  geom_line(data = ipd_curves) +
  ylim(c(0, NA)) + 
  theme_bw() +
  xlab("Age of vaccination") +
  ylab("Incidence (cases per 100,000)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

ggsave(here("output","incidence_plot.pdf"),
       plot = incidence_plot,
       width = 7, height = 5, unit="in")

#plot scaled incidence
ipd_scaled <- ipd %>% group_by(serogroup) %>%
  mutate(p = incidence/sum(incidence))

incidence_plot <- ggplot(data = ipd_scaled,
                         aes(x = agey,
                             y = p,
                             color = serogroup)) + 
  geom_line() +
  theme_bw() +
  xlab("Age (years)") +
  ylab("Scaled incidence") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

ggsave(here("output","scaled_plot.pdf"),
       plot = incidence_plot,
       width = 7, height = 5, unit="in")

#generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- inner_join(ipd_curves, countries_df, by = "agey") %>%
  dplyr::filter(serogroup != "All serotypes") %>%
  mutate(cases = incidence/1e5*ntotal)

#estimate vaccine impact against all IPD serotypes

source(here("script", "metacurve.R"))
#VE_table <- read_csv(here("output","VE_table.csv"))

VE_table <- add_row(VE_table,
                    Study = "None", rate = 0, `Half-life` = Inf)


### make a list check it twice

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

scenarios <- list(`1` = data.frame(Study.waning = "Andrews (2012)",
                                   Study.VE     = "Andrews (2012)"),
                  `2` = data.frame(Study.waning = "Djennad (2018)",
                                   Study.VE     = "Djennad (2018)"),
                  `3` = data.frame(Study.waning = "Andrews (2012)"),
                  `4` = data.frame(Study.waning = "Djennad (2018)")) %>%
  bind_rows(.id = "scenario") %>%
  mutate(age_dep = is.na(Study.VE)) %>%
  crossing(Vac.age = seq(55, 85, by = 5),
           serogroup = c("PPV23",
                         "PCV13")) %>%
  mutate(Study.waning = ifelse(serogroup == "PCV13" & scenario %in% c(2,4),
                               "Andrews (2012)",
                               Study.waning)) %>%
  left_join(dplyr::select(VE_table,
                          Study.waning = Study,
                          rate)) %>%
  left_join(dplyr::select(VE_table,
                          Study.VE = Study,
                          VE)) %>%
  mutate(age_dep = is.na(Study.VE)) %>% 
  mutate(VE = ifelse(age_dep,
                     initial_VE(Vac.age, serogroup, age_dep),
                     VE)) %>%
  mutate(rate  = ifelse(serogroup == "PCV13" & scenario %in% c(1,3),
                        0,
                        rate)) %>%
  mutate(delay = ifelse(serogroup == "PCV13" & scenario %in% c(2,4),
                        5,
                        0))

# we want the curve to be at VE if age >= vac.age + delay
# when delay > 0, we want to subtract delay off

VE_by_Vac.age <- 
  scenarios %>%
  inner_join(Cases) %>%
  mutate(agey_since = agey - Vac.age) %>%
  mutate(agey_since = ifelse(delay == 0,
                             agey_since,
                             pmax(0, agey_since - delay)
  )) %>% 
  mutate(Vaccine_Efficacy = VE*exp(rate*(1 + agey_since))) %>%
  mutate(value = ifelse(agey < Vac.age, 0, Vaccine_Efficacy)) %>%
  mutate(Impact = value*cases) 

impact_by_age_to_plot <- 
  VE_by_Vac.age %>%
  group_by(serogroup,
           Study.waning,
           Study.VE,
           age_dep, delay,
           rate,
           Vac.age, Country,
           scenario) %>%
  summarise(Impact = sum(Impact)) %>%
  mutate(Waning = case_when(rate         == 0 ~ "No waning",
                            Study.waning == "Andrews (2012)" ~ "Fast waning",
                            Study.waning == "Djennad (2018)" ~ "Slow waning"),
         Waning = ifelse(delay > 0,
                         paste(Waning, sprintf("\n(%i years' delay)", delay)),
                         Waning)) 

impact_by_age_to_plot_max <- 
  impact_by_age_to_plot %>%
  group_by_at(.vars = vars(-c(Vac.age, Impact))) %>%
  filter(Impact == max(Impact))

impact_by_age_plot <- 
  ggplot(data = impact_by_age_to_plot,
         aes(x = Vac.age, y= Impact,
             color  = factor(age_dep),
             group = interaction(Waning, age_dep, serogroup, delay,
                                 Country))) +
  geom_line() + 
  facet_grid(Country ~ serogroup + Waning,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (expected total cases averted)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                     palette = "Set1") +
  geom_vline(data = impact_by_age_to_plot_max,
             aes(xintercept = Vac.age,
                 color  = factor(age_dep)),
             lty = 2, alpha = 0.5)

ggsave(filename = "output/impact_by_vac_age.pdf", 
       plot = impact_by_age_plot,
       width = 7, height = 4, units = "in")

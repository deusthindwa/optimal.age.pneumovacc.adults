#load the require packages
if (!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(char = c("tidyverse", "curl", "Hmisc", "here",
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
           nls.control(maxiter=200),
           start = list(log_alpha = log(alpha0),
                        log_beta  = log(beta0),
                        log_theta = log(theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))

ipd_curves <-
  ipd_models %>%
  map_df(~mutate(ipd_x, incidence = predict(object = .x, 
                                            newdata = ipd_x)),
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

#generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- inner_join(ipd_curves, countries_df, by = "agey") %>%
  dplyr::filter(serogroup != "All serotypes") %>%
  mutate(cases = incidence/1e5*ntotal)

#estimate vaccine impact against all IPD serotypes

source(here("script", "metacurve.R"))
#VE_table <- read_csv(here("output","VE_table.csv"))

initial_VE <- function(age, serogroup){
  # age-dependent vaccine efficacy at time of vaccination
  # may be superseded
  dplyr::case_when(serogroup == "PCV13" ~ 0.7,
                   age >= 55 ~ 0.17,
                   age >= 65 ~ 0.25,
                   age >= 75 ~ 0.35,
                   TRUE    ~ NA_real_)
}

VE_by_Vac.age <- 
  dplyr::filter(VE_table,
                Study %in% c("Andrews (2012)",
                             "Djennad (2018)",
                             "All")) %>%
  mutate(Study = fct_recode(factor(Study), c("Pooled" = "All")),
         Study = fct_inorder(Study)) %>%
  crossing(Vac.age = seq(55, 85, by=5)) %>%
  crossing(Cases) %>%
  mutate(VE = initial_VE(Vac.age, serogroup)) %>%
  mutate(Vaccine_Efficacy = VE*exp(rate*(agey - Vac.age + 1))) %>%
  mutate(value = ifelse(agey < Vac.age, 0, Vaccine_Efficacy)) %>%
  mutate(Impact = value*cases) 

impact_by_age_plot <- 
  VE_by_Vac.age %>%
  group_by(Study, serogroup, Vac.age, Country) %>%
  summarise(Impact = sum(Impact)) %>%
  ggplot(data = ., aes(x = Vac.age, y= Impact)) +
  geom_line(aes(color = serogroup)) + 
  facet_grid(Country ~ Study,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (expected total cases averted)") +
  scale_color_brewer(palette= "Set2",
                     name = "Vaccine serotypes"
                     ) +
  theme(legend.position = "bottom")

ggsave(filename = "output/impact_by_vac_age.pdf", 
       plot = impact_by_age_plot,
       width = 7, height = 4, units = "in")

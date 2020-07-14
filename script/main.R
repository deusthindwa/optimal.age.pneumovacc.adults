#load the require packages
require(pacman)
pacman::p_load(char = c("tidyverse", "curl", "Hmisc"))

#load the IPD cases
ipd <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_ipd_incid.csv")) 
pop.ew <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_total_pop.csv")) 
pop.mw <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/MW_total_pop.csv")) 

source("pops.R")

ipd <- mutate(ipd, agey = readr::parse_number(substr(agegroup,1,2)))

#estimate the rest of parameters using a simple linear model
theta0 <- min(ipd$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=ipd)  
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2] 

#Initial parameter values
start <- list(alpha=alpha0, beta=beta0, theta=theta0)

#fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm

IPD_models <- IPD %>% 
  split(.$serogroup) %>%
  map(~nls(data = .x, 
           incidence ~ exp(log_alpha) * exp(exp(log_beta)*agey) + exp(log_theta),
           nls.control(maxiter=200),
           start = list(log_alpha = log(alpha0),
                        log_beta  = log(beta0),
                        log_theta = log(theta0))))

IPD_x <- data.frame(agey = seq(55, 90, by = 1))

IPD_curves <-
  IPD_models %>%
  map_df(~mutate(IPD_x, incidence = predict(object = .x, 
                                            newdata = IPD_x)),
         .id = "serogroup")

#plot fitted curves with backward or forward extrapolation
incidence_plot <- ggplot(data = IPD,
                         aes(x = agey,
                             y = incidence,
                             color = serogroup)) + 
  geom_point(size=2.5) + 
  geom_line(data = IPD_curves) +
  ylim(c(0, NA)) + 
  theme_bw() +
  xlab("Age of vaccination") +
  ylab("Incidence (cases per 100,000)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

ggsave("output/incidence_plot.pdf", plot = incidence_plot,
       width = 7, height = 5, unit="in")

#generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- inner_join(IPD_curves, POP, by = "agey") %>%
  mutate(cases = incidence/1e5*ntotal)

#estimate vaccine impact against all IPD serotypes
Ve      =  0.5
half    =  5
Vac.age = 55

VE_table <- read_csv("output/VE_table.csv")

impact <- filter(VE_table,
                 Study %in% c("Andrews (2012)",
                              "Djennad (2018)")) %>%
  split(.$Study) %>%
  map_df(~mutate(Cases, VE = 0.01*.$VE*exp(.$rate*(agey - Vac.age + 1))),
         .id = "Study") %>%
  mutate(Impact = VE*cases)


VE_by_Vac.age <- 
  filter(VE_table,
       Study %in% c("Andrews (2012)",
                    "Djennad (2018)",
                    "All")) %>%
  mutate(Study = fct_recode(factor(Study), c("Pooled" = "All"))) %>%
  crossing(Vac.age = seq(55, 85, by=5)) %>%
  rowwise %>%
  group_split %>%
  map_df(~mutate(Cases, VE = 0.01*.x$VE*exp(.x$rate*(agey - .x$Vac.age + 1))) %>%
        mutate(Vac.age = .x$Vac.age,
               Study   = .x$Study)) %>%
  mutate(value = ifelse(agey < Vac.age, 0, VE)) %>%
  mutate(Impact = value*cases) 

impact_by_age_plot <- 
  VE_by_Vac.age %>%
  group_by(Study, serogroup, Vac.age) %>%
  summarise(Impact = sum(Impact)) %>%
  ggplot(data = ., aes(x = Vac.age, y= Impact)) +
  geom_line(aes(color = Study)) + 
  facet_grid(. ~ serogroup) +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (expected total cases averted)") +
  scale_color_brewer(palette= "Set2") +
  theme(legend.position = "bottom")

ggsave(filename = "output/impact_by_vac_age.pdf", plot = impact_by_age_plot,
       width = 7, height = 3.5, units = "in")

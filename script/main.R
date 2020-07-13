#load the require packages
require(pacman)
pacman::p_load(char = c("tidyverse", "curl", "Hmisc"))

#load the IPD cases
IPD <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_ipd_incid.csv")) 
POP <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_total_pop.csv")) 

IPD <- mutate(IPD, agey = readr::parse_number(substr(agegroup,1,2)))

#estimate the rest of parameters using a simple linear model
theta0 <- min(IPD$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=IPD)  
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
  ylab("Incidence (%)") +
  theme(legend.position = "bottom") 

ggsave("output/incidence_plot.pdf", plot = incidence_plot,
       width = 7, height = 5, unit="in")

#generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- inner_join(IPD_curves, POP) %>%
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
  map_df(~mutate(Cases, VE = 0.01*.$A*exp(.$B*(agey - Vac.age + 1))),
         .id = "Study") %>%
  mutate(Impact = VE*cases)

ggplot(data = impact,
       aes(x = agey, y = Impact)) +
  geom_line() +
  facet_grid(Study ~ serogroup)

Cases <- Cases %>% mutate(VE=c(rep(0,Vac.age-55),Ve*2^(-1/half*1:(36+55-Vac.age)))) %>% mutate(Impact=VE*all.cases)
plot(Cases$agey, Cases$Impact,type="l",col="blue", ylim=c(0,3e+06),xlim=c(Vac.age,90))
points(Cases$Impact %>% sum)

for(i in c(0.5)){
  for(j in c(60,65,70,75,80,85)){
    Ve=i
    half=5
    Vac.age=j
    Cases <- Cases %>% mutate(VE=c(rep(0,Vac.age-55),Ve*2^(-1/half*1:(36+55-Vac.age)))) %>% mutate(Impact=VE*all.cases)
    lines(Cases$agey, Cases$Impact,col=topo.colors(j,alpha=1),xlim=c(Vac.age,90), lwd=2)
  }
}

Ve =0.2 #first year Ve against VT IPD (currently age independent)
half = 5 #half life of Vein years assuming exponential decay
Vac.age = 60 #year of vaccination

# model
Cases <- Cases %>% 
  mutate(VE=c(rep(0,Vac.age-55),Ve*2^(-1/half*1:(36+55-Vac.age)))) %>%
  mutate(Impact=VE*all.cases)

Cases$Impact %>% sum %>% print
plot(Cases$agey, Cases$Impact)
lines(Cases$agey, Cases$Impact)





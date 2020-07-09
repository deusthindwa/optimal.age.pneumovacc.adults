#load the require package
require(tidyverse)
require(curl)
require(Hmisc)

#load the IPD cases
IPD <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_ipd_cases.csv")) 
POP <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/optimal.age.pneumovacc.adults/master/data/EW_total_pop.csv")) 

#add age increment
IPD$agey <- if_else(IPD$agegroup=="55-59",55,
                      if_else(IPD$agegroup=="60-64",60,
                              if_else(IPD$agegroup=="65-69",65,
                                      if_else(IPD$agegroup=="70-74",70,
                                              if_else(IPD$agegroup=="75-79",75,
                                                      if_else(IPD$agegroup=="80-84",80,85))))))

#estimate the rest of parameters using a simple linear model
theta0 <- min(IPD$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=IPD)  
alpha0 <- exp(coef(model0)[1])
beta0 <- coef(model0)[2] 

#Initial parameter values
start <- list(alpha=alpha0, beta=beta0, theta=theta0)

#fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
model.all <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="All IPD cases")), nls.control(maxiter=200))
model.pcv13 <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="PCV13 IPD")), nls.control(maxiter=200))
model.ppv23 <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="PPV23 IPD")), nls.control(maxiter=200))

#IPD incidence extrapolated for each age-year
Cases <- data_frame(agey=seq(from=55,to=90,by=1)) %>%
                    mutate(all.incid=predict(model.all,list(agey=agey))) %>%
                    mutate(pcv13.incid=predict(model.pcv13,list(agey=agey))) %>%
                    mutate(ppv23.incid=predict(model.ppv23,list(agey=agey)))
Cases <- merge(Cases,POP)
Cases$all.cases <- Cases$all.incid*Cases$ntotal/100000
Cases$pcv13.cases <- Cases$pcv13.incid*Cases$ntotal/100000
Cases$ppv23.cases <- Cases$ppv23.incid*Cases$ntotal/100000

#plot fitted curves with backward or forward extrapolation
ggplot() + 
  geom_point(aes(x=IPD$agey,y=IPD$incidence,color=IPD$serogroup), size=2.5) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.all,list(agey=seq(from=55,to=90,by=1)))), color='#F8766D', size=1) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.pcv13,list(agey=seq(from=55,to=90,by=1)))), color='#00BA38', size=1) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.ppv23,list(agey=seq(from=55,to=90,by=1)))), color='#619CFF', size=1) + 
  ylim(0,125) + xlim(55,90) +
  theme_bw()


# to include baseline VE decay with age

# input
Ve = .5 #first year Ve against VT IPD (currently age independent)
half = 5 #half life of Vein years assuming exponential decay
Vac.age = 65 #year of vaccination

# model
data.frame(age = 55:100,
           Pop = 910000 - cumsum(rep(19,46))*1000,  # pup UK like
           IPD.inc = 1.1^(1:46)) %>%                # IPD cases made up
  mutate(IPD.cases = Pop*IPD.inc) %>%
  mutate(VE = c(rep(0,Vac.age-55),Ve*2^(-1/half*1:(46+55-Vac.age)))) %>%
  mutate(Impact = VE * IPD.cases) -> dt

dt$Impact %>% sum %>% print


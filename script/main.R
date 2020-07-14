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
model.all <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="All IPD cases")), nls.control(maxiter=200))
model.pcv13 <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="PCV13 IPD")), nls.control(maxiter=200))
model.ppv23 <- nls(incidence ~ alpha * exp(beta*agey) + theta, start=start, data=na.omit(subset(IPD,serogroup=="PPV23 IPD")), nls.control(maxiter=200))

#plot fitted curves with backward or forward extrapolation
ggplot() + 
  geom_point(aes(x=IPD$agey,y=IPD$incidence,color=IPD$serogroup), size=2.5) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.all,list(agey=seq(from=55,to=90,by=1)))), color='#F8766D', size=1) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.pcv13,list(agey=seq(from=55,to=90,by=1)))), color='#00BA38', size=1) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.ppv23,list(agey=seq(from=55,to=90,by=1)))), color='#619CFF', size=1) + 
  ylim(0,125) + xlim(55,90) +
  theme_bw()

#generate IPD cases from total pop and IPD incidence annually
Cases <- data_frame(agey=seq(from=55,to=90,by=1)) %>% mutate(all.incid=predict(model.all,list(agey=agey))) %>%
  mutate(pcv13.incid=predict(model.pcv13,list(agey=agey))) %>% mutate(ppv23.incid=predict(model.ppv23,list(agey=agey)))
Cases <- merge(Cases,POP)
Cases$all.cases <- Cases$all.incid*Cases$ntotal
Cases$pcv13.cases <- Cases$pcv13.incid*Cases$ntotal
Cases$ppv23.cases <- Cases$ppv23.incid*Cases$ntotal

#estimate vaccine impact against all IPD serotypes
Ve=0.5
half=5
Vac.age=55
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




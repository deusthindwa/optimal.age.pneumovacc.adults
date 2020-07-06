require(tidyverse)

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
  mutate(VE = c(rep(0,Vac.age-55),Ve*2^(-1/half*1:(46+55-Vac.age)))) 
  mutate(Impact = VE * IPD.cases) -> dt

dt$Impact %>% sum %>% print


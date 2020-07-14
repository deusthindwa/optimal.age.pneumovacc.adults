
#plot observed incidence and fitted curves with backward or forward extrapolation

st.all <- subset(ipd,serogroup=="All serotypes")
st.pcv13 <- subset(ipd,serogroup=="PCV13")
st.ppv23 <- subset(ipd,serogroup=="PPV23")
serogroup <- c("all"="red", "ppv23"="blue", "pcv13"="green")

A <- ggplot() + 
  geom_line(aes(st.all$agey,y=st.all$incidence/sum(st.all$incidence), color='all'), size=0.8) +
  geom_line(aes(st.ppv23$agey,y=st.ppv23$incidence/sum(st.ppv23$incidence), color='ppv23'), size=0.8) + 
  geom_line(aes(st.pcv13$agey,y=st.pcv13$incidence/sum(st.pcv13$incidence), color='pcv13'), size=0.8) + 
  scale_x_continuous(breaks=seq(55,90,5)) + 
  scale_y_continuous(breaks=seq(0,0.5,0.05)) + 
  theme_bw() + 
  labs(title="A", x="Age (years)", y="Observed incidence", color="Serogroup") +
  scale_color_manual(values=serogroup) +
  theme(axis.text.x=element_text(face="bold", size=10, color="black"), axis.text.y=element_text(face="bold", size=10, color="black"))

B <- ggplot() + 
  geom_point(aes(x=ipd$agey,y=ipd$incidence,color=ipd$serogroup), size=2.5) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.all,list(agey=seq(from=55,to=90,by=1)))), color='red', size=0.8) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.ppv23,list(agey=seq(from=55,to=90,by=1)))), color='blue', size=0.8) + 
  geom_line(aes(x=seq(from=55,to=90,by=1),y=predict(model.pcv13,list(agey=seq(from=55,to=90,by=1)))), color='green', size=0.8) + 
  scale_x_continuous(breaks = seq(55,90,5)) + 
  scale_y_continuous(breaks = seq(0,125,20)) + 
  theme_bw() + 
  labs(title="B", x="Age (years)", y="incidence per 100,000 population", color="Serogroup") +
  theme(axis.text.x=element_text(face="bold", size=10, color="black"), axis.text.y=element_text(face="bold", size=10, color="black"))

print(ggarrange(A,B,ncol=2,common.legend=TRUE,legend="right"))
remove(A,B,st.all,st.ppv23,st.pcv13,serogroup)

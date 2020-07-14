
#ggplot comparing % populations in England/Wales versus Malawi

countries <- c("England/Wales"="blue", "Malawi"="red")

ggplot() +
  geom_line(aes(x=pop.ew$agey, y=pop.ew$ntotal*100/sum(pop.ew$ntotal), color="England/Wales")) + 
  geom_line(aes(x=pop.mw$agey, y=pop.mw$ntotal*100/sum(pop.mw$ntotal), color="Malawi")) + 
  scale_x_continuous(breaks = seq(55,90,5)) + 
  scale_y_continuous(breaks = seq(0,6.5,0.5)) + 
  theme_bw() + 
  labs(title="", x="Age (years)", y="Population %", color="Countries") +
  scale_color_manual(values=countries) +
  theme(axis.text.x=element_text(face="bold", size=10, color="black"), axis.text.y=element_text(face="bold", size=10, color="black"))
  
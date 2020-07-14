
#ggplot comparing % populations in England/Wales versus Malawi

countries <- c("England/Wales"="blue", "Malawi"="red")

countries_df <- list(`England/Wales` = pop.ew,
                     `Malawi`        = pop.mw) %>%
  bind_rows(.id = "Country") %>%
  group_by(Country) %>%
  mutate(p = ntotal/sum(ntotal))

countries_plot <- ggplot(data = countries_df,
       aes(x = agey,
           y = p)) +
  geom_col(aes(fill = Country),
           position = position_dodge()) +
  scale_x_continuous(breaks = scales::breaks_width(5)) + 
  scale_y_continuous(labels  = scales::percent,
                     limits = c(0, NA)) + 
  theme_bw() + 
  labs(title="", x="Age (years)", y="Population %", color="Countries") +
  #scale_color_manual(values=countries) +
  theme(axis.text=element_text(face="bold", size=10, color="black"),
        legend.position = "bottom")

ggsave(filename = here("output","countries.pdf"), 
       plot = countries_plot,
       width = 7, height = 3.5, units = "in")

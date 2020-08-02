# make plots

#backward or forward extrapolation incidence plot


incidence_plot <-  
  ggplot(data = ipd_curves, aes(x = agey, y = `50%`,
                                color = serogroup,
                                fill  = serogroup)) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`),
              alpha = 0.2,
              color = NA,) +
  #facet_wrap(~serogroup, nrow = 1) +
  ylim(c(0, NA)) + 
  theme_bw() +
  xlab("Age (years)") +
  ylab("Incidence (cases per 100,000)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_point(data = ipd, aes(y = incidence))



ggsave(here("output","incidence_plot.png"),
       plot = incidence_plot,
       width = 7, height = 5, unit="in", dpi = 300)


#scaled incidence plot
scaled_incidence_plot <-
  ggplot(data = ipd_scaled,
         aes(x = agey,
             y = p,
             color = serogroup)) + 
  geom_line() +
  theme_bw() +
  xlab("Age (years)") +
  ylab("Observed incidence") +
  ylim(c(0, NA)) +
  scale_color_brewer(palette = "Dark2", guide = FALSE)


ggsave(here("output","scaled_plot.png"),
       plot = scaled_incidence_plot,
       width = 7, height = 5, unit="in", dpi = 300)

incidence_plots <- 
  scaled_incidence_plot +
  incidence_plot +
  patchwork::plot_layout(nrow = 1,
                         guides = "collect") &
  theme(legend.position='bottom')

ggsave(here("output", "incidence_plots.png"),
       plot = incidence_plots,
       width = 7, height = 4, unit="in", dpi = 300)



impact_by_age_to_plot_ <-
  impact_by_age_to_plot %>%
  nest(data = c(sim, Impact)) %>%
  mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q)

#vaccine impact plot
impact_by_age_plot <- 
  ggplot(data = impact_by_age_to_plot_,
         aes(x = Vac.age, y= `50%`,
             color  = factor(age_dep),
             group = interaction(Waning, age_dep, serogroup, delay,
                                 Country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = factor(age_dep)),
              color = NA,
              alpha = 0.2) +
  facet_grid(Country ~ serogroup + Waning,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0,NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (expected total cases averted)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                     palette = "Set1") + 
  scale_fill_brewer(name = "Age dependent\nvaccine efficacy",
                    palette = "Set1")

ggsave(filename = "output/impact_by_vac_age.png", 
       plot = impact_by_age_plot,
       width = 7, height = 4, units = "in", dpi = 300)


#impact per 10000 older adults vaccinated
impact_validated_plot <- 
  ggplot(data = impact_validated,
         aes(x = Vac.age, y= `50%`,
             color  = factor(age_dep),
             group = interaction(Waning, age_dep, serogroup, delay,
                                 Country))) +
  geom_line() + 
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = factor(age_dep)),
              color = NA,
              alpha = 0.2) +
  facet_grid(Country ~ serogroup + Waning,
             scales = "free_y") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA)) +
  xlab("Vaccination Age") +
  ylab("Impact (cases averted per 10K older adults vaccinated)") +
  theme(legend.position = "bottom") +
  scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                     palette = "Set1") +
  scale_fill_brewer(name = "Age dependent\nvaccine efficacy",
                    palette = "Set1") +
  theme(panel.grid.minor.y = element_blank())



ggsave(filename = "output/impact_validated_by_vac_age.png", 
       plot = impact_validated_plot,
       width = 7, height = 4, units = "in", dpi = 300)

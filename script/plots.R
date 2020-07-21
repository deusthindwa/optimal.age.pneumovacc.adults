# make plots


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


scaled_incidence_plot <- ggplot(data = ipd_scaled,
                         aes(x = agey,
                             y = p,
                             color = serogroup)) + 
    geom_line() +
    theme_bw() +
    xlab("Age (years)") +
    ylab("Scaled incidence") +
    theme(legend.position = "bottom") +
    ylim(c(0, NA)) +
    scale_color_brewer(palette = "Dark2")

ggsave(here("output","scaled_plot.pdf"),
       plot = scaled_incidence_plot,
       width = 7, height = 5, unit="in")

impact_by_age_plot <- 
    ggplot(data = impact_by_age_to_plot,
           aes(x = Vac.age, y= Impact,
               color  = factor(age_dep),
               group = interaction(Waning, age_dep, serogroup, delay,
                                   Country))) +
    geom_line() + 
    facet_grid(Country ~ serogroup + Waning,
               scales = "free_y") +
    theme_bw() +
    scale_y_continuous(limits = c(0,NA)) +
    xlab("Vaccination Age") +
    ylab("Impact (expected total cases averted)") +
    theme(legend.position = "bottom") +
    scale_color_brewer(name = "Age dependent\nvaccine efficacy",
                       palette = "Set1") +
    geom_vline(data = impact_by_age_to_plot_max,
               aes(xintercept = Vac.age,
                   color  = factor(age_dep)),
               lty = 2, alpha = 0.5)

ggsave(filename = "output/impact_by_vac_age.pdf", 
       plot = impact_by_age_plot,
       width = 7, height = 4, units = "in")

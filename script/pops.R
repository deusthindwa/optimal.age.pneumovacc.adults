
#ggplot comparing % populations in England/Wales versus Malawi

pop.ew <- read_csv(here("data", "EW_total_pop.csv"))
pop.mw <- read_csv(here("data", "MW_total_pop.csv")) 

countries <- c("England/Wales"="blue", "Malawi"="red")

pop.totals <- list(`England/Wales` = 56286961 + 3152879, # mid-2019
                   `Malawi`        = 18628747) %>%
  map_df(.id = "Country", ~data.frame(N = .x))

use.pop.totals <- FALSE

countries_df <- list(`England/Wales` = pop.ew,
                     `Malawi`        = pop.mw) %>%
  bind_rows(.id = "Country") %>%
  group_by(Country) %>%
  inner_join(pop.totals) %>%
  mutate(p = ntotal/(use.pop.totals*N + (1-use.pop.totals)*sum(ntotal)))

countries_plot <- ggplot(data = countries_df,
       aes(x = agey,
           y = p)) +
  geom_col(aes(fill = Country),
           position = position_dodge()) +
  scale_x_continuous(breaks  = seq(50, 90, by = 10),
                     labels  = function(x){gsub(pattern = "90", replacement = "90+", x = x)}) + 
  scale_y_continuous(labels  = scales::percent,
                     limits = c(0, NA)) + 
  theme_bw() + 
  labs(title="", x="Age (years)", 
       y=paste0("Age as share of\n",
               ifelse(use.pop.totals,"total national","55+" ),
               " population"),
       color="Countries") +
  #scale_color_manual(values=countries) +
  theme(axis.text=element_text(face="bold", size=10, color="black"),
        legend.position = "bottom")

ggsave(filename = here("output","countries.png"), 
       plot = countries_plot,
       width = 7, height = 3.5, units = "in", dpi = 300)

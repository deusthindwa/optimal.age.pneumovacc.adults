dat <- list(`Andrews (2012)` = 
              list(`0-2` = c(48, 32, 60),
                   `2-5` = c(21, 03, 06),
                   `5-Inf` = c(15, -3, 30)),
            `Wright (2013)` = 
              list(`0-5` = c(-9, -119, 43),
                   `5-10` = c(38.3, -6, 64),
                   `10-Inf` = c(-21, -137, 35)),
            `Rudnick (2013)` = 
              list(`0-5` = c(41.3, 20, 57),
                   `5-Inf` = c(34.1, 6, 54)),
            `Guttierez (2014)` = 
              list(`0-5` = c(44.5, 19, 62),
                   `5-Inf` = c(32.5, -6, 57)),
            `Djennad (2018)` =
              list(`0-2` = c(41, 23, 54),
                   `2-5` = c(34, 16, 48),
                   `5-Inf` = c(23, 12, 32)))

dat_ <- lapply(X = dat, FUN = function(x){
  map_df(x, ~set_names(.x, c("Mean", "Min", "Max")),
         .id = "Ages")}) %>%
  bind_rows(.id = "Study") %>%
  separate(Ages, into = c("xmin", "xmax")) %>%
  mutate_at(.vars = vars(xmin, xmax),
            .funs = parse_number) %>%
  mutate(xmax = ifelse(is.na(xmax), 20, xmax))

df <- dat_ %>% filter(Study != "Wright (2013)") %>%
  rename(y = "Mean")

f <- function(parms, df){
  
  g <- function(parms, xmin, xmax){
    mean(parms[1]*exp(parms[2]*seq(xmin, xmax, by = 1)))
  }
  
  df <- df %>% rowwise %>% mutate(ynew = g(parms, xmin, xmax)) %>% ungroup
  
  return( sum( (df$y - df$ynew)^2 ))
  
  
}

#ans <- optim(par = c(50, -0.5), fn = f, df = df)

ans_by_study <-
  df %>% split(.$Study) %>%
  map(~optim(par = c(50, -0.5), fn = f, df = .x, hessian = TRUE)) 

simulate_from_ans <- function(x, nsim = 1e4, newdata){
  
  f <- function(x){
    mutate(x, fit = X1*exp(X2*t))
  }
  dat <- data.frame(rmvnorm(n = nsim,
                            mean = x$M,
                            sigma = x$V)) %>%
    mutate(sim = 1:n()) %>%
    crossing(newdata) %>%
    f 
  
  return(dat)
}

df_from_study <- ans_by_study %>%
  map(~list(M = .x$par,
            V = solve(.x$hessian))) %>%
  map_df(~simulate_from_ans(.x, nsim = Nsims, 
                            newdata = data.frame(t = seq(0,20))),
         .id = "Study") %>%
  rename(VE = X1,
         rate = X2) %>%
  mutate(VE = VE/100)

df_by_study_q <- df_from_study %>%
  nest(data = -c(Study, t)) %>%
  mutate(Q = map(data, ~quantile(.x$fit, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data)

ans_by_study_parms_from_model <-
  df_from_study %>%
  distinct(Study, VE, sim, rate) %>%
  mutate(rate = -rate) %>%
  gather(key, value, VE, rate) %>%
  nest(data = c(sim, value)) %>%
  mutate(Q = map(data, ~quantile(.x$value, probs = c(0.025, 0.5, 0.975)))) %>%
  unnest_wider(Q) %>%
  select(-data) %>%
  group_by(Study, key) %>%
  transmute(value = sprintf("%0.2f (%0.2f, %0.2f)", `50%`, `2.5%`, `97.5%`)) %>%
  spread(key, value) %>%
  select(Study, `Initial efficacy` = VE, `Waning rate` = rate) %>%
  write_csv(here("output", "waning_rates.csv"))

ans_by_study_parms <-  # these are estimated parameters from optim
  ans_by_study %>%
  map("par") %>%
  map_df(.id = "Study",
         .f = ~data.frame(A = .x[[1]], B = .x[[2]]))

VE_plot <- ggplot(data=df) +
  geom_segment(aes(x=xmin, xend = xmax,
                   y = y, yend = y)) +
  # stat_function(fun = function(parms, x){parms[1]*exp(parms[2]*x)},
  #               args = list(parms = c(ans$par)), inherit.aes = FALSE,
  #               alpha = 0.5) +
  ylim(c(0, NA)) +
  geom_line(data= df_by_study_q, aes(x=t, y = `50%`)) +
  geom_ribbon(data = df_by_study_q, aes(x = t,
                                         ymin = `2.5%`,
                                         ymax = `97.5%`),
              alpha = 0.25) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  xlab("Years since vaccination (t)") +
  ylab("Vaccine efficacy (VE, %)") +
  scale_color_brewer(palette = "Set1") +
  guides(col = guide_legend(ncol = 2)) +
  facet_wrap(~Study) +
  geom_text(data = ans_by_study_parms,
            x = 15, y = 50, parse = T,
            aes(label = paste("VE == ", round(A,1),
                              "*e^{", round(B,3),
                              "*t}"))) 

ggsave("output/VE_plot.png",
       plot = VE_plot,
       width = 7, height = 7, unit="in", dpi = 300)

VE_table <- ans_by_study_parms %>%
  #add_row(Study = "All", A = A, B = B) %>%
  mutate(`Half-life` = -log(2)/B) %>%
  rename(VE = A, rate = B) %>%
  mutate(VE = VE/100)

write_csv(x = VE_table, path = "output/VE_table.csv")

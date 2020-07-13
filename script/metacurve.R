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
    rename(y = Mean)

f <- function(parms, df){
    
    g <- function(parms, xmin, xmax){
        mean(parms[1]*exp(parms[2]*seq(xmin, xmax, by = 1)))
    }
    
    df <- df %>% rowwise %>% mutate(ynew = g(parms, xmin, xmax)) %>% ungroup
    
    return( sum( (df$y - df$ynew)^2 ))
    
    
}

ans <- optim(par = c(50, -0.5), fn = f, df = df)

A <- round(ans$par[1],2)
B <- round(ans$par[2],4)

df %>%
    #gather(key, value, Young, Old) %>%
    ggplot(data=., aes(x=xmin, xend = xmax,
                       color = Study)) +
    geom_segment(aes(group = Study,
                     y = y, yend = y)) +
    stat_function(fun = function(parms, x){parms[1]*exp(parms[2]*x)},
                  args = list(parms = c(ans$par)), inherit.aes = FALSE) +
    ylim(c(0, NA)) +
    theme_bw() +
    theme(legend.position = 'bottom') +
    xlab("Years since vaccination (t)") +
    ylab("Vaccine efficacy (VE, %)") +
    scale_color_brewer(palette = "Set1") +
    guides(col = guide_legend(ncol = 2)) +
    ggtitle(label =  bquote(VE == .(A)*e^{.(B)*t}))

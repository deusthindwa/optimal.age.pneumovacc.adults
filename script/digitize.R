library(pacman)
p_load(char = c("digitize", "mgcv", "tidyverse"))

# the image has been digitzed
# a <- digitize("data/gr1_lrg.jpg")
# write_csv(a, "data/gr1_lrg.csv")

a <- read_csv("data/gr1_lrg.csv") %>%
    mutate(x = round( (x-9)*16/7, 1),
           y = round(y, 0)/100)

X_spline <- splines::bs(a$x, knots = c(1,4,8))

a_spline <- glm(data=a, y ~ X_spline, family = quasibinomial())

a_spline <- gam(data=a, y ~ s(x, bs="cs", 
                              k = 3,
                              fx = TRUE),
                knots = list(x = c(1,4,8)),
                family = quasibinomial())

bind_cols(a,
          data.frame(
              predict(object = a_spline, 
                      type = "link", se= T))) %>%
    mutate(L = fit - 1.96*se.fit,
           U = fit + 1.96*se.fit) %>%
    mutate_at(.vars = vars(fit, L, U),
              .funs = boot::inv.logit) %>%
    ggplot(data = ., aes(x=x)) +
    geom_ribbon(alpha = 0.5, aes(ymin = L, ymax = U),
                fill = NA, color = "black", lty=2) +
    geom_line(aes(y=fit)) +
    geom_point(aes(y=y)) +
    theme_bw() +
    xlab("Years since vaccination") +
    ylab("Vaccination efficacy") +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA))


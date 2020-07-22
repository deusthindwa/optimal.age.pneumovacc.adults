# estimate the rest of parameters using a simple linear model
theta0 <- min(ipd$incidence,na.rm=TRUE)*0.5  
model0 <- lm(log(incidence-theta0) ~ agey, data=ipd)  
alpha0 <- exp(coef(model0)[1])
beta0  <- coef(model0)[2]


# fit nonlinear (weighted) least-squares estimates of the parameters using Gauss-Newton algorithm
# we parameterise in terms of log-rates to ensure the estimates are positive
ipd_models <- ipd %>% 
    split(.$serogroup) %>%
    purrr::map(~nls(data = .x, 
                    incidence ~ exp(log_alpha) * exp(exp(log_beta)*agey) + exp(log_theta),
                    nls.control(maxiter = 200),
                    start = list(log_alpha = log(alpha0),
                                 log_beta  = log(beta0),
                                 log_theta = log(theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))

ipd_curves <- ipd_models %>%
    purrr::map_df(~dplyr::mutate(ipd_x,
                                 incidence = predict(object = .x,
                                                     newdata = ipd_x)),
                  .id = "serogroup")


# calculate scaled incidence
ipd_scaled <- ipd %>% 
    dplyr::group_by(serogroup) %>%
    dplyr::mutate(p = incidence/sum(incidence))

# generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- dplyr::inner_join(ipd_curves, countries_df, by = "agey") %>%
    dplyr::filter(serogroup != "All serotypes") %>%
    dplyr::mutate(cases = incidence/1e5*ntotal)
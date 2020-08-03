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
                    incidence ~ exp(log_alpha) * exp(beta*agey) + (theta),
                    nls.control(maxiter = 200),
                    start = list(log_alpha = (alpha0),
                                 beta  = (beta0),
                                 theta = (theta0))))

ipd_x <- data.frame(agey = seq(55, 90, by = 1))


simulate_from_model <- function(x, nsim = 1e4, newdata){
    V <- vcov(x)
    M <- coef(x)
    
    #browser()
    
    #f <- eval(rhs(x$call$formula))
    
    dat <- data.frame(rmvnorm(n = nsim, mean = M, sigma = V)) %>%
        mutate(sim = 1:n()) %>%
        crossing(newdata) %>%
        mutate(fit = predict(x, .)) %>%
        dplyr::select(fit, sim, one_of(names(newdata)))
    
    return(dat)
    
}

summarise_from_model <- function(x, probs = c(0.025, 0.5, 0.975)){
    nest(x, data = -agey) %>%
        mutate(Q = map(data, ~quantile(.x$fit,
                                       probs = probs))) %>%
        unnest_wider(Q) %>%
        dplyr::select(-data) %>%
        return
}

ipd_mc <-
    ipd_models %>%
    map(~simulate_from_model(.x, newdata = ipd_x, nsim = Nsims)) 

ipd_curves <- ipd_mc %>%
    map_df(summarise_from_model, .id = "serogroup") %>%
    mutate(serogroup = fct_inorder(factor(serogroup))) 

# calculate scaled incidence
ipd_scaled <- ipd %>% 
    dplyr::group_by(serogroup) %>%
    dplyr::mutate(p = incidence/sum(incidence))

# generate IPD cases from total pop and IPD incidence annually
# table 7
Cases <- dplyr::inner_join(bind_rows(ipd_mc, .id="serogroup"), 
                           countries_df, by = "agey") %>%
    dplyr::filter(serogroup != "All serotypes") %>%
    dplyr::mutate(cases = fit/1e5*ntotal, Vac.age = agey)

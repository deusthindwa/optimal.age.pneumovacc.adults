# OPA waning

opa <- readr::read_csv(here("data","opa_ve.csv")) %>%
    dplyr::filter(serogroup == "PPV23", measure == "OPA GMT")

opa_plot <- ggplot(data=opa, aes(x=age, y = prop)) +
    geom_point() 

opa_theta0 <- min(opa$prop,na.rm=TRUE)*0.5  
opa_model0 <- lm(log(prop-opa_theta0) ~ age, data=opa)  
opa_alpha0 <- (coef(opa_model0)[1])
opa_beta0  <- (coef(opa_model0)[2])

opa_model <- nls(data = opa, 
                 prop ~ alpha * exp(beta*age),
                 nls.control(maxiter=200),
                 start = list(alpha = (opa_alpha0),
                              beta  = -(opa_beta0)))

opa_pred <- data.frame(age = seq(55,95,by=1)) %>%
    mutate(prop = predict(opa_model, newdata = .))

opa_plot + geom_line(data = opa_pred) +
    theme_bw() +
    xlab("Age (years)") +
    ylab("OPA GMT")

OPA_VE <- coef(opa_model)[["beta"]]

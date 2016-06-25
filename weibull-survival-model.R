## ----setup, eval = T, results = 'hide', echo = F-------------------------
knitr::opts_chunk$set(message = FALSE,
                      warning = FALSE,
                      fig.width = 6,
                      fig.height = 4
                      )

## ----load-packages, eval = T, echo = F-----------------------------------
library(httr)
library(readr)
library(cgdsr)
library(purrr)
library(dplyr)
library(assertthat)
library(ggplot2)
library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))
library(shinystan)
library(gridExtra)
library(ggfortify)
library(scales)
library(biostan)

## ----locate-stan-file----------------------------------------------------
if (!require(biostan))
    devtools::install_github('jburos/biostan')
library(biostan)
stan_file <- system.file('stan', 'weibull_survival_null_model.stan', package =  'biostan')

## ----print-stan-code-----------------------------------------------------
biostan::print_stan_file(stan_file)

## ----edit-stan-file, eval = F--------------------------------------------
#  if (interactive())
#      file.edit(stan_file)

## ----view-data-block-----------------------------------------------------
print_stan_file(stan_file, section = 'data')

## ----print-model-block---------------------------------------------------
print_stan_file(stan_file, section = 'model')

## ----print-parameters-block----------------------------------------------
print_stan_file(stan_file, section = 'transformed parameters')

## ----single-value-alpha--------------------------------------------------
alpha_raw <- 0.2
tau_al <- 10
log_alpha <- alpha_raw * tau_al
alpha <- exp(log_alpha)
print(alpha)

## ----dist-alpha----------------------------------------------------------
alpha_raw <- rnorm(1000, 0, 1)
tau_al <- 10
log_alpha <- alpha_raw * tau_al
alpha <- exp(log_alpha)
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha)) + 
    geom_density() + 
    scale_x_log10(labels = scientific)

## ----dist-alpha-vs-raw---------------------------------------------------
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha, y = alpha_raw)) + 
    geom_density2d() + 
    scale_x_log10(labels = scientific)

## ----sim-data-function---------------------------------------------------
sim_data <- function(alpha, mu, Nobs, Ncen) {
    observed_data <- data.frame(os_status = rep_len('DECEASED', Nobs),
                                os_months = rweibull(n = Nobs, alpha, exp(-(mu)/alpha)),
                                stringsAsFactors = F
                                )
    
    censored_data <- data.frame(os_status = rep_len('LIVING', Ncen),
                                os_months = runif(Ncen) * rweibull(Ncen, alpha, exp(-(mu)/alpha)),
                                stringsAsFactors = F
                                )
    
    return(observed_data %>% bind_rows(censored_data))
}

## ----simulated-data------------------------------------------------------
test_alpha <- 0.8
test_mu <- -3

## sample sizes from TCGA blca data
test_nobs <- 179 
test_ncen <- 230

## test these inputs for arbitrary values of alpha & mu
simulated_data <- 
    sim_data(alpha = test_alpha,
                 mu = test_mu,
                 Nobs = test_nobs,
                 Ncen = test_ncen
                 ) 
head(simulated_data)

## ----sim-km-curve--------------------------------------------------------
## plot KM curve from simulated data
simulated_data <- 
    simulated_data %>%
    dplyr::mutate(os_deceased = os_status == 'DECEASED')

autoplot(survival::survfit(Surv(os_months, os_deceased) ~ 1,
                      data = simulated_data
                      ), conf.int = F) + 
    ggtitle('Simulated KM curve')

## ----review-data---------------------------------------------------------
print_stan_file(stan_file, section = 'data')

## ----stan-data-----------------------------------------------------------
observed_data <- simulated_data %>%
    dplyr::filter(os_status == 'DECEASED')

censored_data <- simulated_data %>%
    dplyr::filter(os_status != 'DECEASED')

stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months
)
rm(censored_data)
rm(observed_data)
str(stan_data)

## ----gen-stan-data-function----------------------------------------------
gen_stan_data <- function(data) {
    observed_data <- data %>%
        dplyr::filter(os_status == 'DECEASED')
    
    censored_data <- data %>%
        dplyr::filter(os_status != 'DECEASED')
    
    stan_data <- list(
        Nobs = nrow(observed_data),
        Ncen = nrow(censored_data),
        yobs = observed_data$os_months,
        ycen = censored_data$os_months
    )
}

## ----first-stan-run, warning = TRUE--------------------------------------
recover_simulated <- 
    rstan::stan(stan_file,
                data = stan_data,
                chains = 4,
                iter = 1000,
                seed = 1328025050
                )
print(recover_simulated)

## ------------------------------------------------------------------------
print_stan_file(stan_file, section = 'parameters')

## ----stan-init-values----------------------------------------------------
gen_inits <- function() {
      list(
        alpha_raw = 0.01*rnorm(1),
        mu = rnorm(1)
      )
}

## ----stanfit-with-inits--------------------------------------------------
recover_simulated2 <- 
    rstan::stan(stan_file,
                data = stan_data,
                chains = 4,
                iter = 1000,
                init = gen_inits
                )
print(recover_simulated2)

## ----traceplot-lp--------------------------------------------------------
rstan::traceplot(recover_simulated2, 'lp__')

## ------------------------------------------------------------------------
rstan::traceplot(recover_simulated2, c('alpha','mu'), ncol = 1)

## ----launch-shinystan, eval = F------------------------------------------
#  if (interactive())
#      shinystan::launch_shinystan(recover_simulated2)

## ----stanfit-only-observed-----------------------------------------------
recover_simulated_obs <- 
    rstan::stan(stan_file,
                data = gen_stan_data(
                    simulated_data %>% dplyr::filter(os_status == 'DECEASED')
                    ),
                chains = 4,
                iter = 1000,
                init = gen_inits
                )
print(recover_simulated_obs)

## ------------------------------------------------------------------------
recover_simulated_cen <- 
    rstan::stan(stan_file,
                data = gen_stan_data(
                    simulated_data %>% dplyr::filter(os_status != 'DECEASED')
                    ),
                chains = 4,
                iter = 1000,
                init = gen_inits
                )
print(recover_simulated_cen)

## ----sim-extract-alpha---------------------------------------------------
pp_alpha <- rstan::extract(recover_simulated2,'alpha')$alpha
pp_mu <- rstan::extract(recover_simulated2,'mu')$mu

## ----sim-post-predict----------------------------------------------------
pp_newdata <- 
    purrr::map2(.x = pp_alpha,
                .y = pp_mu,
                .f = ~ sim_data(alpha = .x, 
                                mu = .y,
                                Nobs = test_nobs,
                                Ncen = test_ncen
                                )
                )

## ----sim-plot-time-to-event----------------------------------------------
ggplot(pp_newdata %>%
           dplyr::bind_rows() %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           bind_rows(simulated_data %>% dplyr::mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~type, ncol = 1)

## ----sim-pp-survdata-----------------------------------------------------
## cumulative survival rate for each posterior draw
pp_survdata <-
    pp_newdata %>%
    purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
    purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
    purrr::map(fortify)

## summarize cum survival for each unit time (month), summarized at 95% confidence interval
pp_survdata_agg <- 
    pp_survdata %>%
    purrr::map(~ dplyr::mutate(., time_group = floor(time))) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarize(surv_mean = mean(surv)
                     , surv_p50 = median(surv)
                     , surv_lower = quantile(surv, probs = 0.025)
                     , surv_upper = quantile(surv, probs = 0.975)
                     ) %>%
    dplyr::ungroup()

## ----sim-plot-ppcheck----------------------------------------------------
## km-curve for test data 
test_data_kmcurve <- 
    fortify(
        survival::survfit(
            Surv(os_months, os_deceased) ~ 1, 
            data = simulated_data %>% 
                dplyr::mutate(os_deceased = os_status == 'DECEASED')
            )) %>%
    dplyr::mutate(lower = surv, upper = surv)

ggplot(pp_survdata_agg %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           dplyr::rename(surv = surv_p50, lower = surv_lower, upper = surv_upper, time = time_group) %>%
           bind_rows(test_data_kmcurve %>% dplyr::mutate(type = 'actual data')),
       aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    xlim(c(0, 200))

## ----pp_predict-function-------------------------------------------------
pp_predict_surv <- function(pp_alpha, pp_mu, Nobs, Ncen,
                            level = 0.9,
                            plot = F, data = NULL,
                            sim_data_fun = sim_data
                            ) {
    pp_newdata <- 
        purrr::map2(.x = pp_alpha,
                    .y = pp_mu,
                    .f = ~ sim_data_fun(alpha = .x, mu = .y,
                                    Nobs = Nobs, Ncen = Ncen
                                    )
                    )
    
    pp_survdata <-
        pp_newdata %>%
        purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
        purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
        purrr::map(fortify)
    
    ## compute quantiles given level 
    lower_p <- 0 + ((1 - level)/2)
    upper_p <- 1 - ((1 - level)/2)
    
    pp_survdata_agg <- 
        pp_survdata %>%
        purrr::map(~ dplyr::mutate(.,
                                   time_group = floor(time))) %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(time_group) %>%
        dplyr::summarize(surv_mean = mean(surv)
                         , surv_p50 = median(surv)
                         , surv_lower = quantile(surv,
                                                 probs = lower_p)
                         , surv_upper = quantile(surv,
                                                 probs = upper_p)
                         ) %>%
        dplyr::ungroup()
    
    if (plot == FALSE) {
        return(pp_survdata_agg)
    } 
    
    ggplot_data <- pp_survdata_agg %>%
           dplyr::mutate(type = 'posterior predicted values') %>%
           dplyr::rename(surv = surv_p50,
                         lower = surv_lower,
                         upper = surv_upper, time = time_group)
    
    if (!is.null(data))
        ggplot_data <- 
            ggplot_data %>% 
            bind_rows(
                fortify(
                    survival::survfit(
                        Surv(os_months, os_deceased) ~ 1, 
                        data = data %>% 
                            dplyr::mutate(
                                os_deceased = os_status == 'DECEASED')
                        )) %>%
                dplyr::mutate(lower = surv,
                              upper = surv, type = 'actual data')
                )
    
    pl <- ggplot(ggplot_data,
                 aes(x = time, group = type, linetype = type)) + 
        geom_line(aes(y = surv, colour = type)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
        
    pl 
}

## ----ppcheck-obs-only----------------------------------------------------
pp_alpha_obs <- extract(recover_simulated_obs, 'alpha')$alpha
pp_mu_obs <- extract(recover_simulated_obs, 'mu')$mu
pl <- pp_predict_surv(pp_alpha = pp_alpha_obs,
                pp_mu = pp_mu_obs,
                Nobs = test_nobs,
                Ncen = 0,
                plot = T,
                data = simulated_data %>% 
                    dplyr::filter(os_status == 'DECEASED')
                )
pl + 
    ggtitle('PP check against simulated data\n(model fit using observed events only)') +
    xlim(NA, 200)

## ----alt-sim-data-function-----------------------------------------------
alt_sim_data <- function(alpha, mu, Nobs, Ncen) {
    
    data <- data.frame(surv_months = rweibull(n = Nobs + Ncen, alpha, exp(-(mu)/alpha)),
                       censor_months = rexp(n = Nobs + Ncen, rate = 1/100),
                       stringsAsFactors = F
                       ) %>%
        dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                          'DECEASED', 'LIVING'
                                          ),
                       os_months = ifelse(surv_months < censor_months,
                                          surv_months, censor_months
                                          )
                       )

    return(data)
}

## ----alt-sim-data--------------------------------------------------------

alt_simulated_data <- alt_sim_data(
    alpha = test_alpha,
    mu = test_mu,
    Ncen = test_ncen,
    Nobs = test_nobs
    )
autoplot(survival::survfit(Surv(os_months, I(os_status == 'DECEASED')) ~ 1,
                           data = alt_simulated_data
                           ))

## ------------------------------------------------------------------------
table(alt_simulated_data$os_status)

## ----recover-alt-sim-----------------------------------------------------
recover_alt_simulated <- rstan::stan(
    file = stan_file,
    data = gen_stan_data(alt_simulated_data),
    chains = 4,
    iter = 1000,
    init = gen_inits
)
print(recover_alt_simulated)

## ------------------------------------------------------------------------
pp_alpha <- rstan::extract(recover_alt_simulated, 'alpha')$alpha
pp_mu <- rstan::extract(recover_alt_simulated, 'mu')$mu
pl <- pp_predict_surv(pp_alpha = pp_alpha,
                pp_mu = pp_mu,
                sim_data_fun = alt_sim_data,
                Nobs = test_nobs, Ncen = test_ncen,
                plot = T, data = alt_simulated_data
                )
pl + ggtitle('KM curve for actual & posterior predicted data\nAfter modifying the simulate-data function')


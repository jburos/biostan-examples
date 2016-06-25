## ----setup, eval = T, results = 'hide', echo = F-------------------------
knitr::opts_chunk$set(message = FALSE,
                      warning = FALSE,
                      fig.width = 6,
                      fig.height = 6
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

## ----example-load-data, eval = FALSE-------------------------------------
#  url <- 'http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id=blca_tcga_all'
#  req <- httr::GET(url)
#  clinical_data <-
#      httr::content(req,
#                    type = 'text/tab-separated-values',
#                    col_names = T,
#                    col_types = NULL
#                    )
#  str(clinical_data)

## ----actual-load-data----------------------------------------------------
mycgds = cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
selected_case_list = 'blca_tcga_all'
clinical_data = cgdsr::getClinicalData(mycgds, selected_case_list)

## ----inspect-data--------------------------------------------------------
str(clinical_data,  no.list = T, vec.len = 2)

## ----initprep-data-------------------------------------------------------
## names to lower case
names(clinical_data) <- tolower(names(clinical_data))

## convert empty strings -> NA values
convert_blank_to_na <- function(x) {
    if (!purrr::is_character(x)) {
        warning('input vector is not character - returning original input')
        return(x)
    } else {
        ifelse(x == '', NA, x)
    }
}
clin_data <- clinical_data %>%
    dplyr::mutate_each(funs = funs(convert_blank_to_na), everything())

## inspect resulting data frame
str(clin_data, vec.len = 2, list.len = 10)

## ----inspect-os-status---------------------------------------------------
clinical_data %>%
    dplyr::filter(is.na(os_status) | os_status == '') %>%
    dplyr::select(os_status, os_months) %>%
    str()

## ----inspect-os-months---------------------------------------------------
clinical_data %>%
    dplyr::filter(!is.na(os_status) & os_status != '') %>%
    dplyr::filter(os_months < 0 | is.na(os_months)) %>%
    dplyr::select(os_status, os_months) %>%
    head()

## ----create-clin-data----------------------------------------------------
clin_data <- 
    clin_data %>%
    dplyr::filter(!is.na(os_status) & os_status != '') %>%
    dplyr::filter(os_months >= 0 & !is.na(os_months))

## confirm 4 fewer observations than original
assert_that(nrow(clin_data) == nrow(clinical_data) - 4)

## ----plot-os-months------------------------------------------------------
ggplot(clin_data,
       aes(x = os_months,
           group = os_status,
           colour = os_status,
           fill = os_status
           )) + 
    geom_density(alpha = 0.5)

## ----km-survival-curve---------------------------------------------------
mle.surv <- 
    survfit(
        Surv(os_months,os_deceased) ~ 1,
        data = clin_data %>%
            dplyr::mutate(os_deceased = os_status == 'DECEASED')
    )
autoplot(mle.surv, conf.int = F) +
    ggtitle('KM survival curve for BLCA cohort')

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
weibull_sim_data <- function(alpha, mu, n) {
    
    data <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu)/alpha)),
                       censor_months = rexp(n = n, rate = 1/100),
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

## ----simulated-data------------------------------------------------------
test_alpha <- 0.8
test_mu <- -3

## sample size from TCGA blca data
test_n <- nrow(clin_data)

## test these inputs for arbitrary values of alpha & mu
simulated_data <- 
    weibull_sim_data(alpha = test_alpha,
                 mu = test_mu,
                 n = test_n
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
                data = gen_stan_data(simulated_data),
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

## ----sim-stanfit-with-inits----------------------------------------------
recover_simulated2 <- 
    rstan::stan(stan_file,
                data = gen_stan_data(simulated_data),
                chains = 4,
                iter = 1000,
                init = gen_inits
                )
print(recover_simulated2)

## ----sim-traceplot-lp, fig.height=3--------------------------------------
rstan::traceplot(recover_simulated2, 'lp__')

## ----sim-traceplot-params------------------------------------------------
rstan::traceplot(recover_simulated2, c('alpha','mu'), ncol = 1)

## ----sim-launch-shinystan, eval = F--------------------------------------
#  if (interactive())
#      shinystan::launch_shinystan(recover_simulated2)

## ----sim-extract-alpha---------------------------------------------------
pp_alpha <- rstan::extract(recover_simulated2,'alpha')$alpha
pp_mu <- rstan::extract(recover_simulated2,'mu')$mu

## ----plot-alpha-vs-test--------------------------------------------------
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density(aes(x = alpha)) + 
    geom_vline(aes(xintercept = test_alpha), colour = 'red') +
    ggtitle('Posterior distribution of alpha\nshowing true value in red')

## ----plot-mu-vs-test-----------------------------------------------------
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density(aes(x = mu)) + 
    geom_vline(aes(xintercept = test_mu), colour = 'red') +
    ggtitle('Posterior distribution of mu\nshowing true value in red')

## ----plot-mu-vs-alpha----------------------------------------------------
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
    geom_density2d(aes(x = alpha, y = mu)) +
    geom_point(aes(x = test_alpha, y = test_mu), colour = 'red', size = 2) +
    ggtitle('Posterior distributions of mu and alpha\nshowing true parameter values in red')

## ------------------------------------------------------------------------
mean(pp_alpha >= test_alpha)

## ------------------------------------------------------------------------
mean(pp_mu >= test_mu)

## ------------------------------------------------------------------------
mean(pp_mu >= test_mu & pp_alpha >= test_alpha)

## ----sim-post-predict----------------------------------------------------
pp_newdata <- 
    purrr::map2(.x = pp_alpha,
                .y = pp_mu,
                .f = ~ weibull_sim_data(alpha = .x, 
                                mu = .y,
                                n = test_n
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
pp_predict_surv <- function(pp_alpha, pp_mu, n,
                            level = 0.9,
                            plot = F, data = NULL,
                            sim_data_fun = weibull_sim_data
                            ) {
    pp_newdata <- 
        purrr::map2(.x = pp_alpha,
                    .y = pp_mu,
                    .f = ~ sim_data_fun(alpha = .x, mu = .y, n = n)
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

## ----fit-model-tgca------------------------------------------------------
wei_fit <- rstan::stan(file = stan_file,
                       data = gen_stan_data(clin_data),
                       iter = 1000,
                       chains = 4,
                       init = gen_inits
                       )
print(wei_fit)

## ----lp-traceplot--------------------------------------------------------
rstan::traceplot(wei_fit, c('lp__'), ncol = 1)

## ----param-traceplot-----------------------------------------------------
rstan::traceplot(wei_fit, c('alpha','mu'), ncol = 1)

## ----launch-shinystan, eval = F------------------------------------------
#  if (interactive())
#      launch_shinystan(wei_fit)

## ----wei-ppchecks--------------------------------------------------------
pl <- pp_predict_surv(pp_alpha = extract(wei_fit,'alpha')$alpha,
                pp_mu = extract(wei_fit,'mu')$mu,
                n = nrow(clin_data),
                data = clin_data,
                plot = T
                ) 
pl + 
    xlim(NA, 150) +
    ggtitle('Posterior predictive checks for NULL weibull model\nfit to TCGA data; showing 90% CI')

## ----summarize-coverage--------------------------------------------------
## summarize 90% CI of predicted event rate for each interval
pp_agg <- pp_predict_surv(pp_alpha = extract(wei_fit,'alpha')$alpha,
                pp_mu = extract(wei_fit,'mu')$mu,
                n = nrow(clin_data)
                )


## summarize observed data into same time_groups
act_agg <- 
    survival::survfit(Surv(os_months, I(os_status == 'DECEASED')) ~ 1,
                             data = clin_data
                             ) %>%
    fortify() %>%
    dplyr::mutate(time_group = floor(time)) %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarise(observed_surv = mean(surv)) %>%
    dplyr::ungroup()

## compute proportion of observed values within 90% ci
act_agg %>%
    dplyr::inner_join(pp_agg, by = 'time_group') %>%
    dplyr::mutate(within_interval = ifelse(observed_surv >= surv_lower & observed_surv <= surv_upper,
                                           1, 0),
                  time_set = cut(time_group, breaks = c(0, 100))
                  ) %>%
    dplyr::group_by(time_set) %>%
    dplyr::summarize(mean(within_interval))


library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(bayestestR)
library(multidplyr)

## Read in the data from Henne & O'Neill (Experiment 1)
d <- read_csv('https://osf.io/download/3psjc/') %>%
    mutate(Cause=Cause/100+.5)

## Plot the data for reference:
## Note that we will use the ordbetareg package to account for
## the prevalence of 1's and 0's in the data.
ggplot(d, aes(x=Cause, fill=Condition)) +
    geom_histogram(binwidth=.01) +
    theme_classic()

## Make a fake dataset of n_participants participants per condition
pred_data <- function(d, n_participants=1) {
    d %>%
        distinct(Condition) %>%
        expand_grid(id=1:n_participants) %>%
        mutate(Cause=NA)
}

## Convert the data d to be used by stan
## d: the dataframe to embed as a list for stan
## n_participants: the number of participants per condition to simulate in X_pred
## prior_only: if true, simulate from the prior (not the posterior)
standata <- function(d, n_participants=1, prior_only=FALSE) {
    mm <- model.matrix(Cause ~ Condition, data=d)
    mm.pred <- model.matrix(Cause ~ Condition, pred_data(d, n_participants))
    list(prior_only=prior_only, N=dim(mm)[1], K=dim(mm)[2],
         X=mm, y=d$Cause,
         N_pred=dim(mm.pred)[1], X_pred=mm.pred)
}

## Simulate n_datasets datasets of causal judgments with n_participants participants per condition
sim_data <- function(n_participants=1, n_datasets=1) {
    m_pred$generate_quantities(fit, data=standata(d, n_participants=n_participants),
                               parallel_chains=4) %>%
        spread_draws(y_pred[.row], ndraws=n_datasets) %>%
        mutate(Condition=pred_data(d, n_participants=n_participants)$Condition[.row]) %>%
        rename(Cause=y_pred, sim=.draw) %>%
        group_by(sim) %>%
        select(-.row, -.chain, -.iteration) %>%
        nest() %>%
        pull(data)
}

## Calculate the BF for the parameter par in model given the prior
log10_bf <- function(model, prior, par) {
    bf_pointnull(as.numeric(model$draws(par, format='draws_matrix')),
                 prior=as.numeric(prior$draws(par, format='draws_matrix')))$log_BF / log(10)
}





## Fit the model in stan
m <- cmdstan_model('ordbeta.stan')
m_pred <- cmdstan_model('ordbeta_pred.stan')
fit <- m$sample(standata(d), parallel_chains=4)
prior <- m$sample(standata(d, prior_only=TRUE), parallel_chains=4)

## Extract BFs
log10_bf(fit, prior, 'b[2]')
log10_bf(fit, prior, 'b_phi[2]')

## Inspect the model fit
sim_data(1000)[[1]] %>%
    ggplot(aes(x=Cause, fill=Condition)) +
    geom_histogram(aes(x=Cause, y=after_stat(density/100)), data=d, binwidth=.01) +
    geom_freqpoly(aes(linetype=Condition, y=after_stat(density/100)), binwidth=.01) +
    xlab('Causal Judgment') + ylab('Probability') +
    theme_classic()






## Set up a cluster to run simulations in parallel
cluster <- new_cluster(6)
cluster_library(cluster, c('tidyverse', 'tidybayes', 'bayestestR'))
cluster_copy(cluster, c('d', 'm', 'm_pred', 'fit', 'prior', 'standata', 'pred_data', 'sim_data'))


## Simulate 1000 datasets of varying size, fit the model, and calculate significance
if (file.exists('power_simulations.rds')) {
    sims <- readRDS('power_simulations.rds')
} else {
    sims <- expand_grid(sim=1:1000,
                        N=c(25, 50, 100)) %>%
        group_by(N) %>%
        mutate(data=sim_data(N[1], n_datasets=n())) %>%
        group_by(N, sim) %>%
        partition(cluster) %>%
        mutate(
            ## fit a model to each dataset
            model=map(data, ~ m$sample(standata(.), refresh=0, show_messages=FALSE)),

            ## check for significant mean difference
            b=map_dbl(model, ~ .$summary()$median[5]),
            b.lower=map_dbl(model, ~ .$summary()$q5[5]),
            b.upper=map_dbl(model, ~ .$summary()$q95[5]),
            log10_BF=map_dbl(model, ~ log10_BF(., prior, 'b[2]')),
            significant=b.lower > 0 & log10_BF > 1,
            
            ## check for significant variance difference
            b_phi=map_dbl(model, ~ .$summary()$median[7]),
            b_phi.lower=map_dbl(model, ~ .$summary()$q5[7]),
            b_phi.upper=map_dbl(model, ~ .$summary()$q95[7]),
            log10_BF_phi=map_dbl(model, ~ log10_BF(., prior, 'b_phi[2]')),
            significant_phi=b_phi.lower > 0 & log10_BF_phi > 1) %>%
    collect()
    saveRDS(sims, 'power_simulations.rds')
}


## Calculate power for each dataset size
sims %>% group_by(N) %>% summarize(power=mean(significant))
sims %>% group_by(N) %>% summarize(power=mean(significant_phi))

## Plot the lower confidence bound of the coefficient of interest
ggplot(sims) +
    geom_histogram(aes(x=b.lower), alpha=.5, fill='orange', binwidth=.05) +
    geom_histogram(aes(x=b_phi.lower), alpha=.5, fill='blue', binwidth=.05) +
    geom_vline(xintercept=0, linetype='dashed') +
    facet_wrap(~ N) +
    theme_classic()
ggsave('power_lower_bound.png', width=6, height=3)

## Plot the BFs for the coefficient of interest
ggplot(sims) +
    geom_histogram(aes(x=log10_BF), alpha=.5, fill='orange', binwidth=.5) +
    geom_histogram(aes(x=log10_BF_phi), alpha=.5, fill='blue', binwidth=.5) +
    geom_vline(xintercept=1, linetype='dashed') +
    facet_wrap(~ N) +
    coord_cartesian(xlim=c(-5, 20)) +
    theme_classic()
ggsave('power_BF.png', width=6, height=3)

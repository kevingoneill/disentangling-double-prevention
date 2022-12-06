library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(rstantools)
library(loo)

## sample from a thresholded normal distribution with bounds [lower, upper]
r_thresh_norm <- function(n, mean=0, sd=1, lower=0, upper=1) {
    pmin(upper, pmax(lower, rnorm(n, mean, sd)))
}

standata <- function(d, prior_only=FALSE) {
    list(N=nrow(d),
         V=length(unique(d$vignette)),
         F=length(unique(d$factor)),
         vignette=d$vignette,
         factor=d$factor,
         cause=d$cause,
         prior_only=prior_only)
}

## Data simulation parameters
params <- expand_grid(factor=factor(c('PC', 'DP'), levels=c('PC', 'DP')),
                      vignette=factor(c('reversal', 'adversarial'),
                                      levels=c('reversal', 'adversarial'))) %>%
    mutate(p.DP=.2, p.PC=.8,
           p.PP=ifelse(vignette=='reversal', .9, .1),
           p.E=p.PC * (1-p.PP*(1-p.DP)),
           ns=ifelse(factor=='PC',
           (1-p.PC) + p.PC*((1-p.PP)*(1-p.DP) + (1-p.PP)*p.DP + p.PP*p.DP),
           (1-p.DP) + p.DP*p.PC),
           ces=ifelse(factor=='PC', sqrt((p.PC*(1-p.PC))/(p.E*(1-p.E))) * (1-p.PP + p.PP*p.DP),
                      sqrt((p.DP*(1-p.DP))/(p.E*(1-p.E))) * p.PC * p.PP),
           sd=.33)
params_long <- params %>%
    select(-factor) %>%
    pivot_longer(p.DP:p.PP, names_to=c('parameter', 'factor'), names_sep='\\.') %>%
    mutate(factor=factor(factor, levels=c('PC', 'DP', 'PP')))


## Compile models
ns <- cmdstan_model('necessity_sufficiency.stan')
ces <- cmdstan_model('counterfactual_effect_size.stan')



## Simulate data from NS model
d_ns <- params %>%
    expand_grid(participant=1:100) %>%
    mutate(cause=r_thresh_norm(n(), ns, sd)) %>%
    select(vignette, factor, cause)

ggplot(d_ns, aes(x=vignette, y=cause, fill=factor)) +
    geom_violin() +
    theme_classic()
ggsave('plots/data_ns.png', width=6, height=4)


## Test NS model on data simulated from itself
fit.ns.ns <- ns$sample(standata(d_ns), parallel_chains=4)
draws.ns.ns <- fit.ns.ns %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d_ns$vignette)[vignette], levels=levels(d_ns$vignette)),
           factor=factor(factor, levels=c(levels(d_ns$factor), 'PP')))

draws.ns.ns %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye() +
    geom_point(aes(y=ns, color=factor),
               position=position_dodge(.1, preserve='single'), data=params) +
    xlab('Vignette') + ylab('Predicted Causal Judgment') +
    ggtitle('NS model fit to data simulated from NS model') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/fit_ns_ns.png', width=6, height=6)

draws.ns.ns %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(scale=.75, side='left', slab_alpha=.33) +
    geom_point(aes(x=vignette, color=factor, y=value), data=params_long,
               position=position_dodge(-.25, preserve='single')) +
    xlab('Vignette') + ylab('Predicted Probability') +
    ylim(0, 1) +
    ggtitle('NS model fit to data simulated from NS model') +
    theme_classic()
ggsave('plots/fit_ns_ns_prob.png', width=6, height=6)



## Test CES model on data simulated from NS model
fit.ces.ns <- ces$sample(standata(d_ns), parallel_chains=4)
draws.ces.ns <- fit.ces.ns %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d_ns$vignette)[vignette], levels=levels(d_ns$vignette)),
           factor=factor(factor, levels=c(levels(d_ns$factor), 'PP')))

draws.ces.ns %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye() +
    geom_point(aes(y=ns, color=factor),
               position=position_dodge(.1, preserve='single'), data=params) +
    xlab('Vignette') + ylab('Predicted Causal Judgment') +
    ggtitle('CES model fit to data simulated from NS model') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/fit_ces_ns.png', width=6, height=6)

draws.ces.ns %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(scale=.75, side='left', slab_alpha=.33) +
    geom_point(aes(x=vignette, color=factor, y=value), data=params_long,
               position=position_dodge(-.25, preserve='single')) +
    xlab('Vignette') + ylab('Predicted Probability') +
    ggtitle('CES model fit to data simulated from NS model') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/fit_ces_ns_prob.png', width=6, height=6)




## Simulate data from CES model
d_ces <- params %>%
    expand_grid(participant=1:100) %>%
    mutate(cause=r_thresh_norm(n(), ces, sd)) %>%
    select(vignette, factor, cause)

ggplot(d_ces, aes(x=vignette, y=cause, fill=factor)) +
    geom_violin() +
    theme_classic()
ggsave('plots/data_ces.png', width=6, height=4)


## Test NS model on data simulated from CES model
fit.ns.ces <- ns$sample(standata(d_ces), parallel_chains=4)
draws.ns.ces <- fit.ns.ces %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d_ns$vignette)[vignette], levels=levels(d_ns$vignette)),
           factor=factor(factor, levels=c(levels(d_ns$factor), 'PP')))

draws.ns.ces %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye() +
    geom_point(aes(y=ces, color=factor),
               position=position_dodge(.1, preserve='single'), data=params) +
    xlab('Vignette') + ylab('Predicted Causal Judgment') +
    ggtitle('NS model fit to data simulated from CES model') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/fit_ns_ces.png', width=6, height=6)

draws.ns.ces %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(scale=.75, side='left', slab_alpha=.33) +
    geom_point(aes(x=vignette, color=factor, y=value), data=params_long,
               position=position_dodge(-.25, preserve='single')) +
    xlab('Vignette') + ylab('Predicted Probability') +
    ggtitle('NS model fit to data simulated from CES model') +
    theme_classic()
ggsave('plots/fit_ns_ces_prob.png', width=6, height=6)


## Test CES model on data simulated from itself
fit.ces.ces <- ces$sample(standata(d_ces), parallel_chains=4)
draws.ces.ces <- fit.ces.ces %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d_ns$vignette)[vignette], levels=levels(d_ns$vignette)),
           factor=factor(factor, levels=c(levels(d_ns$factor), 'PP')))

draws.ces.ces %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye() +
    geom_point(aes(y=ces, color=factor),
               position=position_dodge(.1, preserve='single'), data=params) +
    xlab('Vignette') + ylab('Predicted Causal Judgment') +
    ggtitle('CES model fit to data simulated from CES model') +
    theme_classic()
ggsave('plots/fit_ces_ces.png', width=6, height=6)

draws.ces.ces %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(scale=.75, side='left', slab_alpha=.33) +
    geom_point(aes(x=vignette, color=factor, y=value), data=params_long,
               position=position_dodge(-.25, preserve='single')) +
    xlab('Vignette') + ylab('Predicted Probability') +
    ggtitle('CES model fit to data simulated from CES model') +
    theme_classic()
ggsave('plots/fit_ces_ces_prob.png', width=6, height=6)




## Compare models using R squared
R2.ns.ns <- bayes_R2(fit.ns.ns$draws('cause_hat', format='draws_matrix'), d_ns$cause)
R2.ces.ns <- bayes_R2(fit.ces.ns$draws('cause_hat', format='draws_matrix'), d_ns$cause)
R2.ns.ces <- bayes_R2(fit.ns.ces$draws('cause_hat', format='draws_matrix'), d_ces$cause)
R2.ces.ces <- bayes_R2(fit.ces.ces$draws('cause_hat', format='draws_matrix'), d_ces$cause)

tibble(R2.ns.ns, R2.ces.ns, R2.ns.ces, R2.ces.ces) %>%
    pivot_longer(R2.ns.ns:R2.ces.ces, names_to=c('measure', 'model', 'data'), names_sep='\\.') %>%
    ggplot(aes(x=value, y=model)) +
    xlab('R squared') +
    stat_halfeye() +
    facet_grid(data ~ ., labeller=label_both) +
    theme_classic()
ggsave('plots/r_squared.png', width=6, height=4)

tibble(R2.ns.ns/R2.ces.ns, R2.ns.ces/R2.ces.ces) %>%
    rename(ns=`R2.ns.ns/R2.ces.ns`,
           ces=`R2.ns.ces/R2.ces.ces`) %>%
    pivot_longer(ns:ces, names_to='data') %>%
    ggplot(aes(x=value, y=data)) +
    xlab('R squared ratio (NS / CES)') +
    stat_halfeye() +
    geom_vline(xintercept=1, linetype='dashed') +
    theme_classic()
ggsave('plots/r_squared_ratio.png', width=6, height=4)


## Compare models using leave-one-out cross validation of log likelihood
loo_compare(list(ns=fit.ns.ns$loo(), ces=fit.ces.ns$loo()))
loo_compare(list(ns=fit.ns.ces$loo(), ces=fit.ces.ces$loo()))

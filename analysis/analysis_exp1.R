library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(bayestestR)
library(rstantools)
library(loo)
library(ggdist)
library(bayesplot)

PALETTE <- c("#E61C38", "#C0C0C0", "#0B7189")
PALETTE2 <- colorblind_pal()(3)[c(3, 2)] ##c("#006494", "#65743A")


## Convert the data d to be used by stan
## d: the dataframe to embed as a list for stan
## prior_only: if true, simulate from the prior (not the posterior)
standata <- function(d, prior_only=FALSE) {
    mm <- model.matrix(cause ~ vignette*factor, data=d)
    mm.pred <- d %>% distinct(vignette, factor) %>% model.matrix(~ vignette*factor, data=.)
    
    list(prior_only=prior_only, N=dim(mm)[1], K=dim(mm)[2],
         X=mm, y=d$cause,
         N_pred=dim(mm.pred)[1], X_pred=mm.pred)
}

## read in the data, exclude based on attention check, and relevel factors
d <- read_csv('data/exp1_processed.csv') %>%
    filter(attention_check=='Yes.') %>%
    mutate(factor=factor(factor, levels=c('productive_factor', 'double_preventer')),
           vignette=factor(vignette, levels=c('base', 'reversal', 'adversarial')),)
nrow(d)   ## number of participants after exclusion (out of 630)


## compile and fit the stan model
m <- cmdstan_model('ordbeta.stan')
fit <- m$sample(standata(d), chains=10, parallel_chains=10, iter_sampling=5000)
prior <- m$sample(standata(d, prior=TRUE), chains=10, parallel_chains=10, iter_sampling=5000)

## check for convergence and print summary
fit$cmdstan_diagnose()
fit$summary(c('cutpoints', 'b', 'b_phi'))
d %>% distinct(vignette, factor) %>% model.matrix(~ vignette*factor, data=.) %>% colnames

## plot posterior predictive distribution
ppc_ecdf_overlay(d$cause, fit$draws('y_pred', format='draws_matrix')[1:100, ])

## extract means per condition
draws <- fit %>%
    gather_draws(e_pred[.row], mu[.row], phi[.row]) %>%
    left_join(d %>% distinct(vignette, factor) %>% mutate(.row=row_number())) %>%
    left_join(prior %>%
              gather_draws(e_pred[.row], mu[.row], phi[.row]) %>%
              rename(.prior=.value)) %>%
    group_by(.variable, vignette)

## plot computed means against raw data
draws %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=vignette, fill=factor)) +
    stat_slab(aes(y=cause, side=ifelse(factor=='productive_factor', 'left', 'right')),
              data=d, position=position_dodge(.2)) +
    stat_pointinterval(aes(y=.value, group=factor, justification=ifelse(factor=='productive_factor', 1, 0)),
                       point_interval=median_hdi, .width=.95, position=position_dodge(.2)) +
    scale_x_discrete(name='', labels=c('Base', 'Reversal', 'Adversarial')) +
    scale_y_continuous(expand=c(0, 0), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_fill_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    ylab('Causal Judgment') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_means.png', width=7, height=5)


## contrast mean judgments of each factor within vignettes
contr.factor <- draws %>%
    compare_levels(.value, by='factor', comparison=list(c('productive_factor', 'double_preventer'))) %>%
    left_join(draws %>% compare_levels(.prior, by='factor', comparison=list(c('productive_factor', 'double_preventer')))) %>%
    mutate(log10_BF=bf_pointnull(.value, .prior)$log_BF / log(10))
contr.factor %>%
    group_by(vignette, .variable, log10_BF) %>%
    median_hdi(.value)


contr.factor %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=vignette, y=.value, fill=vignette)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95, show.legend=FALSE) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Base', 'Reversal', 'Adversarial')) +
    scale_fill_manual(values=c('grey80', PALETTE2)) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_factor_epred.png', width=6, height=5)

contr.factor %>%
    filter(.variable=='mu') %>%
    ggplot(aes(x=vignette, y=.value, fill=vignette)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95, show.legend=FALSE) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Base', 'Reversal', 'Adversarial')) +
    scale_fill_manual(values=c('grey80', PALETTE2)) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_factor_mu.png', width=6, height=5)

contr.factor %>%
    filter(.variable=='phi') %>%
    ggplot(aes(x=vignette, y=.value, fill=vignette)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95, show.legend=FALSE) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Base', 'Reversal', 'Adversarial')) +
    scale_fill_manual(values=c('grey80', PALETTE2)) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_factor_phi.png', width=6, height=5)


contr.interaction <- contr.factor %>%
    group_by(.variable) %>%
    compare_levels(.value, by='vignette', comparison='control') %>%
    left_join(contr.factor %>% group_by(.variable) %>% compare_levels(.prior, by='vignette', comparison='control')) %>%
    mutate(log10_BF=bf_pointnull(.value, .prior)$log_BF / log(10),
           vignette=factor(vignette, levels=c('reversal - base', 'adversarial - base'), labels=c('reversal', 'adversarial')))
contr.interaction %>%
    group_by(vignette, .variable, log10_BF) %>%
    median_hdi(.value)


contr.interaction %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Reversal', 'Adversarial')) +
    ylab('Double Prevention Effect Compared to Base Case\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_interaction_epred.png', width=6, height=5)

contr.interaction %>%
    filter(.variable=='mu') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Reversal', 'Adversarial')) +
    ylab('Double Prevention Effect Compared to Base Case\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_interaction_mu.png', width=6, height=5)

contr.interaction %>%
    filter(.variable=='phi') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Reversal', 'Adversarial')) +
    ylab('Double Prevention Effect Compared to Base Case\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment1_contrast_interaction_phi.png', width=6, height=5)





## fit NS & CES model
ns <- cmdstan_model('necessity_sufficiency.stan')
ces <- cmdstan_model('counterfactual_effect_size.stan')

d.models <- d %>% filter(vignette != 'base') %>% mutate(vignette=factor(vignette, levels=c('reversal', 'adversarial')))
d.stan <- list(N=nrow(d.models), V=length(unique(d.models$vignette)), F=length(unique(d.models$factor)),
               vignette=d.models$vignette, factor=d.models$factor, cause=d.models$cause, prior_only=FALSE)

fit.ns <- ns$sample(d.stan, chains=10, parallel_chains=10, iter_sampling=5000)
fit.ces <- ces$sample(d.stan, chains=10, parallel_chains=10, iter_sampling=5000)

fit.ns
fit.ces


draws.ns <- fit.ns %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d.models$vignette)[vignette], levels=levels(d.models$vignette)),
           factor=factor(factor, levels=c('PC', 'DP', 'PP')))
draws.ces <- fit.ces %>%
    spread_draws(p_PC, p_DP, p_PP[vignette], K_PC[vignette], K_DP[vignette]) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(vignette=factor(levels(d.models$vignette)[vignette], levels=levels(d.models$vignette)),
           factor=factor(factor, levels=c('PC', 'DP', 'PP')))

## Plot model fit
draws.ns %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value, color=factor),
                           ##justification=ifelse(factor=='productive_factor', 1, 0)),
                       point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15), show.legend=FALSE,
                       data=draws %>% filter(vignette!='base', .variable=='e_pred') %>%
                           mutate(factor=factor(factor, labels=c('PC', 'DP')))) +
    scale_x_discrete(name='', labels=c('Reversal', 'Adversarial')) +
    ylab('Predicted Causal Judgment') +
    scale_fill_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    scale_color_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    ggtitle('NS Model Fit') +
    ylim(0, 1) +
    theme_classic(base_size=18)
ggsave('plots/experiment1_ns.png', width=6, height=6)

draws.ces %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=vignette, y=value, fill=factor)) +
    stat_halfeye(point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value, color=factor,
                           justification=ifelse(factor=='productive_factor', 1, 0)),
                       point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15), show.legend=FALSE,
                       data=draws %>% filter(vignette!='base', .variable=='e_pred') %>% mutate(factor=factor(factor, labels=c('PC', 'DP')))) +
    scale_x_discrete(name='', labels=c('Reversal', 'Adversarial')) +
    ylab('Predicted Causal Judgment') +
    scale_fill_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    scale_color_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    ggtitle('CES Model Fit') +
    ylim(0, 1) +
    theme_classic(base_size=18)
ggsave('plots/experiment1_ces.png', width=6, height=6)


## Plot model-estimated sampling probabilities
draws.ns %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=factor, y=value, fill=vignette, side=vignette)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95,
                 position=position_dodge(.2), show.legend=c(size=FALSE, side=FALSE)) +
    scale_x_discrete(name='', labels=c('Productive\nFactor', 'Double\nPreventer', 'Possible\nPreventer')) +
    scale_side_mirrored(start='bottomleft') +
    ylab('Predicted Probability of Imagining Event\n(Necessity-Sufficiency Model)') +
    scale_fill_manual(name='', values=PALETTE2, labels=c('Reversal', 'Adversarial')) +
    scale_y_continuous(limits=c(0, 1), expand=c(0,0), labels=c('0', '.25', '.5', '.75', '1')) +
    theme_classic(base_size=18)
ggsave('plots/experiment1_ns_prob.png', width=8, height=6)


draws.ces %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=factor, y=value, fill=vignette, side=vignette)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95,
                 position=position_dodge(.2), show.legend=c(size=FALSE, side=FALSE)) +
    scale_x_discrete(name='', labels=c('Productive\nFactor', 'Double\nPreventer', 'Possible\nPreventer')) +
    scale_side_mirrored(start='bottomleft') +
    ylab('Predicted Probability of Imagining Event\n(Counterfactual Effect Size Model)') +
    scale_fill_manual(name='', values=PALETTE2, labels=c('Reversal', 'Adversarial')) +
    scale_y_continuous(limits=c(0, 1), expand=c(0,0), labels=c('0', '.25', '.5', '.75', '1')) +
    theme_classic(base_size=18)
ggsave('plots/experiment1_ces_prob.png', width=8, height=6)



## compare models using LOO-IC
loo_compare(list(ns=fit.ns$loo(), ces=fit.ces$loo()))

## Compare models using R-squared
R2.ns <- bayes_R2(fit.ns$draws('cause_hat', format='draws_matrix'), d.models$cause)
R2.ces <- bayes_R2(fit.ces$draws('cause_hat', format='draws_matrix'), d.models$cause)

tibble(NS=R2.ns, CES=R2.ces) %>%
    pivot_longer(NS:CES, names_to='model') %>%
    ggplot(aes(x=value, y=model)) +
    xlab('R squared') + ylab('') + xlim(0, NA) +
    stat_halfeye() +
    theme_classic()
ggsave('plots/experiment1_r2.png', width=6, height=4)

tibble(ratio=R2.ns/R2.ces) %>%
    ggplot(aes(x=ratio)) +
    xlab('R squared Ratio (NS / CES)') + ylab('') +
    stat_halfeye() +
    geom_vline(xintercept=1, linetype='dashed') +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank())
ggsave('plots/experiment1_r2_ratio.png', width=6, height=4)



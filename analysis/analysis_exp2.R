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
    mm <- model.matrix(cause ~ condition*factor, data=d)
    mm.pred <- d %>% distinct(condition, factor) %>% model.matrix(~ condition*factor, data=.)
    
    list(prior_only=prior_only, N=dim(mm)[1], K=dim(mm)[2],
         X=mm, y=d$cause,
         V=length(levels(d$vignette)), K_V=dim(mm)[2], v=d$vignette, X_V=mm,
         N_pred=dim(mm.pred)[1], X_pred=mm.pred, X_V_pred=mm.pred)
}


d <- read_csv('data/exp2_processed.csv') %>%
    filter(attention_check=='Yes.') %>%
    mutate(factor=factor(factor, levels=c('productive_factor', 'double_preventer')),
           condition=factor(condition, levels=c('reversal', 'adversarial')),
           vignette=factor(vignette))
nrow(d)   ## number of participants after exclusion (out of 2100)


m <- cmdstan_model('ordbeta_vignette.stan')

if (file.exists('experiment2.rds')) {
    fit <- readRDS('experiment2.rds')
} else {
    fit <- m$sample(standata(d), chains=10, parallel_chains=10, iter_sampling=5000, adapt_delta=.95)
    fit$save_object('experiment2.rds')
}
prior <- m$sample(standata(d, prior=TRUE), chains=10, parallel_chains=10, iter_sampling=5000)

## check for convergence and print summary
fit$cmdstan_diagnose()
fit$summary(c('cutpoints', 'b', 'b_phi'))
d %>% distinct(condition, factor) %>% model.matrix(~ condition*factor, data=.) %>% colnames


fit$summary(c('sd_v', 'Omega')) %>% print(n=100)  # vignette-level effects

## plot posterior predictive distribution
ppc_ecdf_overlay(d$cause, fit$draws('y_pred', format='draws_matrix')[1:100, ])


## extract means per condition
draws <- fit %>%
    gather_draws(e_pred[.row], mu_logit[.row], phi_log[.row]) %>%
    left_join(d %>% distinct(condition, factor) %>% mutate(.row=row_number())) %>%
    left_join(prior %>%
              gather_draws(e_pred[.row], mu_logit[.row], phi_log[.row]) %>%
              rename(.prior=.value)) %>%
    group_by(.variable, condition, factor)


## plot computed means against raw data
draws %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=condition, fill=factor)) +
    stat_slab(aes(y=cause, side=ifelse(factor=='productive_factor', 'left', 'right')),
              data=d, position=position_dodge(.5)) +
    stat_pointinterval(aes(y=.value, group=factor,
                           justification=ifelse(factor=='productive_factor', 1, 0)),
                       point_interval=median_hdi, .width=.95, position=position_dodge(.5)) +
    scale_x_discrete(name='', labels=c('Reversal', 'Adversarial')) +
    scale_y_continuous(expand=c(0, 0), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_fill_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    ylab('Causal Judgment') +
    theme_classic(base_size=18)
ggsave('plots/experiment2_means.png', width=6, height=5)


## contrast mean judgments of each factor within condition
contr.factor <- draws %>%
    compare_levels(.value, by='factor',
                   comparison=list(c('productive_factor', 'double_preventer'))) %>%
    left_join(draws %>%
              compare_levels(.prior, by='factor',
                             comparison=list(c('productive_factor', 'double_preventer')))) %>%
    mutate(BF=exp(bf_pointnull(.value, .prior)$log_BF),
           P=as.numeric(pd_to_p(pd(.value))))

contr.factor %>%
    group_by(condition, .variable, BF, P) %>%
    median_hdi(.value)


contr.factor %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    xlim(-.5, .5) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18) +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
ggsave('plots/experiment2_contrast_factor_epred.png', width=5, height=5)

contr.factor %>%
    filter(.variable=='mu_logit') %>%
    ggplot(aes(y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    xlim(-.5, .5) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18) +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
ggsave('plots/experiment2_contrast_factor_mu.png', width=5, height=5)

contr.factor %>%
    filter(.variable=='phi_log') %>%
    ggplot(aes(y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    xlim(-.5, .5) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic() +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
ggsave('plots/experiment2_contrast_factor_phi.png', width=5, height=5)




contr.interaction <- contr.factor %>%
    group_by(.variable) %>%
    compare_levels(.value, by='condition', comparison='control') %>%
    left_join(contr.factor %>% group_by(.variable) %>%
              compare_levels(.prior, by='condition', comparison='control')) %>%
    mutate(BF=exp(bf_pointnull(.value, .prior)$log_BF),
           P=as.numeric(pd_to_p(pd(.value))))

contr.interaction %>%
    group_by(condition, .variable, BF, P) %>%
    median_hdi(.value)

contr.interaction %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=condition, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Adversarial - Reversal')) +
    ylab('Difference in Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment2_contrast_interaction_epred.png', width=6, height=5)

contr.interaction %>%
    filter(.variable=='mu_logit') %>%
    ggplot(aes(x=condition, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Adversarial - Reversal')) +
    ylab('Difference in Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment2_contrast_interaction_mu.png', width=6, height=5)

contr.interaction %>%
    filter(.variable=='phi_log') %>%
    ggplot(aes(x=condition, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Condition', labels=c('Adversarial - Reversal')) +
    ylab('Difference in Double Prevention Effect on Variance\n(Productive Factor - Double Preventer)') +
    theme_classic(base_size=18)
ggsave('plots/experiment2_contrast_interaction_phi.png', width=6, height=5)


## Vignette-level effects
draws.vignette <- fit %>%
    gather_draws(e_pred_v[.row, vignette], mu_v_logit[.row, vignette], phi_v_log[.row, vignette]) %>%
    left_join(d %>% distinct(condition, factor) %>% mutate(.row=row_number())) %>%
    left_join(prior %>%
              gather_draws(e_pred_v[.row, vignette], mu_v_logit[.row, vignette],
                           phi_v_log[.row, vignette]) %>%
              rename(.prior=.value)) %>%
    mutate(vignette=levels(d$vignette)[vignette],
           .variable=case_when(.variable=='e_pred_v' ~ 'e_pred',
                               .variable=='mu_v_logit' ~ 'mu_logit',
                               .variable=='phi_v_log' ~ 'phi_log')) %>%
    group_by(.variable, condition, factor, vignette)

## plot computed means against raw data
draws.vignette %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=condition, fill=factor)) +
    stat_slab(aes(y=cause, side=ifelse(factor=='productive_factor', 'left', 'right')),
              data=d, position=position_dodge(.5)) +
    stat_pointinterval(aes(y=.value, group=factor, justification=ifelse(factor=='productive_factor', 1, 0)),
                       point_interval=median_hdi, .width=.95, position=position_dodge(.5)) +
    scale_x_discrete(name='', labels=c('Reversal', 'Adversarial')) +
    scale_fill_manual(name='', labels=c('Productive\nFactor', 'Double\nPreventer'), values=PALETTE) +
    facet_grid(~ vignette, labeller=as_labeller(str_to_title)) +
    ylab('Causal Judgment') +
    theme_classic(base_size=18) +
    theme(axis.text.x=element_text(angle=-30, hjust=0))
ggsave('plots/experiment2_means_vignette.png', width=12, height=5)


contr.factor.vignette <- draws.vignette %>%
    compare_levels(.value, by='factor',
                   comparison=list(c('productive_factor', 'double_preventer'))) %>%
    left_join(draws.vignette %>%
              compare_levels(.prior, by='factor',
                             comparison=list(c('productive_factor', 'double_preventer')))) %>%
    mutate(BF=exp(bf_pointnull(.value, .prior)$log_BF),
           P=as.numeric(pd_to_p(pd(.value))))

contr.factor.vignette %>%
    filter(.variable=='mu_logit') %>%
    group_by(condition, vignette, .variable, BF, P) %>%
    median_hdi(.value) %>%
    print(n=30)


contr.factor.vignette %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=vignette, y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_factor_epred_vignette.png', width=10, height=5)

contr.factor.vignette %>%
    filter(.variable=='mu_logit') %>%
    ggplot(aes(x=vignette, y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_factor_mu_vignette.png', width=10, height=5)


contr.factor.vignette %>%
    filter(.variable=='phi_log') %>%
    ggplot(aes(x=vignette, y=.value, side=condition,
               fill=condition)) +
    stat_halfeye(point_interval=median_hdi, .width=.95,
                 position=position_dodge(.1)) +
    geom_hline(yintercept=0, linetype='dashed') +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    scale_fill_manual(values=PALETTE2, name='', labels=str_to_title) +
    scale_side_mirrored(start='bottomleft') + guides(side=FALSE) +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_factor_phi_vignette.png', width=12, height=5)




contr.interaction.vignette <- contr.factor.vignette %>%
    group_by(.variable, vignette) %>%
    compare_levels(.value, by='condition', comparison='control') %>%
    left_join(contr.factor.vignette %>% group_by(.variable, vignette) %>%
              compare_levels(.prior, by='condition', comparison='control')) %>%
    mutate(BF=exp(bf_pointnull(.value, .prior)$log_BF),
           P=as.numeric(pd_to_p(pd(.value))))

contr.interaction.vignette %>%
    group_by(vignette, condition, .variable, BF, P) %>%
    median_hdi(.value)

contr.interaction.vignette %>%
    filter(.variable=='e_pred') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    coord_cartesian(ylim=c(-.5, .5)) +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    ylab('Difference in Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_interaction_epred_vignette.png', width=6, height=5)

contr.interaction.vignette %>%
    filter(.variable=='mu_logit') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    coord_cartesian(ylim=c(-3, 3)) +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    ylab('Difference in Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_interaction_mu_vignette.png', width=6, height=5)

contr.interaction.vignette %>%
    filter(.variable=='phi_log') %>%
    ggplot(aes(x=vignette, y=.value)) +
    ##stat_slab(aes(y=.prior)) +
    stat_halfeye(point_interval=median_hdi, fill='skyblue', .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    coord_cartesian(ylim=c(-3, 3)) +
    scale_x_discrete(name='Vignette', labels=str_to_title) +
    ylab('Difference in Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    theme_classic()
ggsave('plots/experiment2_contrast_interaction_phi_vignette.png', width=6, height=5)








## fit vignette-level models with shift/scale
ns <- cmdstan_model('necessity_sufficiency_vignette.stan')
ces <- cmdstan_model('counterfactual_effect_size_vignette.stan')

d.stan <- list(N=nrow(d), C=length(unique(d$condition)),
               V=length(unique(d$vignette)), F=length(unique(d$factor)),
               condition=d$condition, vignette=d$vignette,
               factor=d$factor, cause=d$cause, prior_only=FALSE)

fit.ns <- ns$sample(d.stan, chains=10, parallel_chains=10, iter_sampling=5000)
fit.ces <- ces$sample(d.stan, chains=10, parallel_chains=10, iter_sampling=5000)

fit.ns$cmdstan_diagnose()
fit.ces$cmdstan_diagnose()


fit.ns
fit.ces


draws.ns <- fit.ns %>%
    spread_draws(p_PC[vignette], p_DP[vignette], p_PP[vignette, condition],
                 K_PC[vignette, condition], K_DP[vignette, condition], shift, scale) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(value=ifelse(parameter=='K', value*scale + shift, value),
           condition=factor(levels(d$condition)[condition], levels=levels(d$condition)),
           vignette=factor(levels(d$vignette)[vignette], levels=levels(d$vignette)),
           factor=factor(factor, levels=c('PC', 'DP', 'PP')))
draws.ces <- fit.ces %>%
    spread_draws(p_PC[vignette], p_DP[vignette], p_PP[vignette, condition],
                 K_PC[vignette, condition], K_DP[vignette, condition], shift, scale) %>%
    pivot_longer(c(starts_with('K_'), starts_with('p_')),
                 names_to=c('parameter', 'factor'), names_sep='_') %>%
    mutate(value=ifelse(parameter=='K', value*scale + shift, value),
           condition=factor(levels(d$condition)[condition], levels=levels(d$condition)),
           vignette=factor(levels(d$vignette)[vignette], levels=levels(d$vignette)),
           factor=factor(factor, levels=c('PC', 'DP', 'PP')))


draws.ns %>%
    filter(parameter=='K') %>%
    pivot_wider(names_from=factor) %>%
    mutate(value=PC-DP) %>%
    ggplot(aes(x=vignette, y=value)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value), point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15),
                       data=contr.factor.vignette %>% filter(.variable=='e_pred')) +
    geom_hline(yintercept=0, linetype='dashed') +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    facet_grid(~ condition, labeller=as_labeller(str_to_title)) +
    theme_classic()
ggsave('plots/experiment2_ns_contrast.png', width=12, height=6)

draws.ces %>%
    filter(parameter=='K') %>%
    pivot_wider(names_from=factor) %>%
    mutate(value=PC-DP) %>%
    ggplot(aes(x=vignette, y=value)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value), point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15),
                       data=contr.factor.vignette %>% filter(.variable=='e_pred')) +
    geom_hline(yintercept=0, linetype='dashed') +
    ylab('Double Prevention Effect\n(Productive Factor - Double Preventer)') +
    facet_grid(~ condition, labeller=as_labeller(str_to_title)) +
    theme_classic()
ggsave('plots/experiment2_ces_contrast.png', width=12, height=6)

draws.ns %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=condition, y=value, fill=factor)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value, color=factor),
                       point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15),
                       show.legend=FALSE,
                       data=draws.vignette %>% filter(.variable=='e_pred') %>%
                           mutate(factor=factor(factor,
                                                levels=c('productive_factor', 'double_preventer'),
                                                labels=c('PC', 'DP')))) +
    scale_x_discrete(name='Condition', labels=str_to_title) +
    ylab('Predicted Causal Judgment') +
    scale_fill_discrete(name='', labels=c('Productive\nFactor', 'Double\nPreventer')) +
    facet_grid(~ vignette, labeller=as_labeller(str_to_title)) +
    ggtitle('NS Model Fit') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/experiment2_ns_vignette2.png', width=12, height=6)

draws.ces %>%
    filter(parameter=='K') %>%
    ggplot(aes(x=condition, y=value, fill=factor)) +
    stat_halfeye(scale=.75, point_interval=median_hdi, .width=.95) +
    stat_pointinterval(aes(y=.value, color=factor),
                       point_interval=median_hdi, .width=.95,
                       position=position_nudge(x=-.15),
                       show.legend=FALSE,
                       data=draws.vignette %>% filter(.variable=='e_pred') %>%
                           mutate(factor=factor(factor,
                                                levels=c('productive_factor', 'double_preventer'),
                                                labels=c('PC', 'DP')))) +
    scale_x_discrete(name='Condition', labels=str_to_title) +
    ylab('Predicted Causal Judgment') +
    scale_fill_discrete(name='', labels=c('Productive\nFactor', 'Double\nPreventer')) +
    facet_grid(~ vignette, labeller=as_labeller(str_to_title)) +
    ggtitle('CES Model Fit') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/experiment2_ces_vignette2.png', width=12, height=6)



## Correlating model predictions with causal judgments
draws.vignette %>%
    filter(.variable=='e_pred') %>%
    left_join(draws.ns %>% filter(parameter=='K', factor!='PP') %>% rename(NS=value) %>%
              mutate(factor=factor(`factor`, labels=levels(draws.vignette$factor)))) %>%
    median_hdi(.value, NS) %>%
    ggplot(aes(x=.value, y=NS)) +
    geom_abline(linetype='dashed', slope=1, intercept=0) +
    geom_errorbarh(aes(xmin=.value.lower, xmax=.value.upper, color=factor)) +
    geom_errorbar(aes(ymin=NS.lower, ymax=NS.upper, color=factor)) +
    geom_point(aes(color=factor)) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE) +
    scale_color_manual(values=PALETTE, name='', labels=c('Productive\nFactor', 'Double\nPreventer')) +
    xlab('Mean Causal Judgment') + ylab('Predicted Causal Judgment\n(NS Model)') +
    theme_bw()
ggsave('plots/experiment2_ns_correlation.png', width=6, height=5)

draws.vignette %>%
    filter(.variable=='e_pred') %>%
    left_join(draws.ces %>% filter(parameter=='K', factor!='PP') %>% rename(CES=value) %>%
              mutate(factor=factor(`factor`, labels=levels(draws.vignette$factor)))) %>%
    median_hdi(.value, CES) %>%
    ggplot(aes(x=.value, y=CES)) +
    geom_abline(linetype='dashed', slope=1, intercept=0) +
    geom_errorbarh(aes(xmin=.value.lower, xmax=.value.upper, color=factor)) +
    geom_errorbar(aes(ymin=CES.lower, ymax=CES.upper, color=factor)) +
    geom_point(aes(color=factor)) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE) +
    scale_color_manual(values=PALETTE, name='', labels=c('Productive\nFactor', 'Double\nPreventer')) +
    xlab('Mean Causal Judgment') + ylab('Predicted Causal Judgment\n(CES Model)') +
    theme_bw()
ggsave('plots/experiment2_ces_correlation.png', width=6, height=5)


## Correlating model probabilities with DP effect
contr.factor.vignette


## Correlating model probability difference with interaction effect
contr.interaction.vignette %>%
    filter(.variable=='e_pred') %>%
    left_join(draws.ns %>%
              filter(parameter=='p', factor=='PP') %>%
              pivot_wider(names_from=condition) %>%
              mutate(prob=reversal-adversarial)) %>%
    median_hdi(.value, prob) %>%
    ggplot(aes(x=.value, y=prob, color=vignette)) +
    geom_errorbarh(aes(xmin=.value.lower, xmax=.value.upper), height=.1) +
    geom_errorbar(aes(ymin=prob.lower, ymax=prob.upper)) +
    geom_point() +
    coord_fixed(xlim=c(-.5, 1), ylim=c(-.5, 1), expand=FALSE) +
    scale_color_discrete(name='Vignette', labels=str_to_title) +
    xlab('Interaction Effect on Causal Judgments') +
    ylab('Sampling Probability Difference (NS Model)') +
    theme_bw()
ggsave('plots/experiment2_ns_interaction.png', width=6, height=5)


contr.interaction.vignette %>%
    filter(.variable=='e_pred') %>%
    left_join(draws.ces %>%
              filter(parameter=='p', factor=='PP') %>%
              pivot_wider(names_from=condition) %>%
              mutate(prob=reversal-adversarial)) %>%
    median_hdi(.value, prob) %>%
    ggplot(aes(x=.value, y=prob, color=vignette)) +
    geom_errorbarh(aes(xmin=.value.lower, xmax=.value.upper)) +
    geom_errorbar(aes(ymin=prob.lower, ymax=prob.upper)) +
    geom_point() +
    coord_fixed(xlim=c(-.1, .65), ylim=c(-.1, .65), expand=FALSE) +
    scale_color_discrete(name='Vignette', labels=str_to_title) +
    xlab('Interaction Effect on Causal Judgments') +
    ylab('Sampling Probability Difference (CES Model)') +
    theme_bw()
ggsave('plots/experiment2_ces_interaction.png', width=6, height=5)




draws.ns %>%
    filter(parameter=='p', vignette=='heartworm', condition=='reversal') %>%
    pivot_wider(names_from='factor') %>%
    mutate(DP_effect=(1 - PC*PP + PC*PP*DP) - (1 - DP + DP*PC),
           PP=ntile(PP, n=5),
           scale=ntile(scale, n=5)) %>%
    ggplot(aes(x=PC, y=DP)) +
    geom_point(aes(color=DP_effect), alpha=.01) +
    facet_grid(scale ~ PP, labeller=label_both) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
    scale_color_continuous(limits=c(-1/3, 0)) +
    theme_classic()


ggsave('NS.png', width=8, height=4)



## Plot model-estimated sampling probabilities
draws.ns %>%
    filter(parameter=='p', factor=='PP') %>%
    pivot_wider(names_from=condition) %>%
    mutate(value=reversal-adversarial) %>%
    ggplot(aes(x=vignette, y=value)) +
    stat_halfeye(scale=.25, side='left', point_interval=median_hdi, .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    ylab('Predicted Probability Difference') +
    theme_classic()
ggsave('plots/experiment2_ns_manip.png', width=12, height=6)

draws.ces %>%
    filter(parameter=='p', factor=='PP') %>%
    pivot_wider(names_from=condition) %>%
    mutate(value=reversal-adversarial) %>%
    ggplot(aes(x=vignette, y=value)) +
    stat_halfeye(scale=.25, side='left', point_interval=median_hdi, .width=.95) +
    geom_hline(yintercept=0, linetype='dashed') +
    ylab('Predicted Probability Difference') +
    theme_classic()
ggsave('plots/experiment2_ces_manip.png', width=12, height=6)




draws.ns %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=factor, y=value, fill=condition)) +
    stat_halfeye(aes(side=condition), show.legend=FALSE,
                 point_interval=median_hdci, .width=.95,
                 scale=.125, slab_alpha=.75, normalize='none',
                 position=position_dodge(.2)) +
    ylab('Predicted Sampling Probability (NS Model)') +
    facet_grid(~ vignette, labeller=as_labeller(str_to_title)) +
    scale_fill_manual(values=PALETTE2, labels=c('Reversal', 'Adversarial')) +
    scale_side_mirrored(start='bottomleft') +
    scale_x_discrete(name='', labels=c('Productive\nCause', 'Double\nPreventer',
                                             'Possible\nPreventer')) +
    ylim(0, 1) +
    theme_classic() +
    theme(axis.title.x=element_blank())

ggsave('plots/experiment2_ns_prob_vignette3.png', width=12, height=6)

draws.ces %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=factor, y=value, fill=condition)) +
    stat_halfeye(aes(side=condition), show.legend=FALSE,
                 point_interval=median_hdci, .width=.95,
                 scale=.1, slab_alpha=.75, normalize='none',
                 position=position_dodge(.2)) +
    ylab('Predicted Sampling Probability (CES Model)') +
    facet_grid(~ vignette, labeller=as_labeller(str_to_title)) +
    scale_fill_manual(values=PALETTE2, labels=c('Reversal', 'Adversarial')) +
    scale_side_mirrored(start='bottomleft') +
    scale_x_discrete(name='', labels=c('Productive\nCause', 'Double\nPreventer',
                                             'Possible\nPreventer')) +
    ylim(0, 1) +
    theme_classic() +
    theme(axis.title.x=element_blank())

draws.ces %>%
    filter(parameter=='p') %>%
    ggplot(aes(x=condition, y=value, fill=factor)) +
    stat_slab(scale=.75, side='left', slab_alpha=.33) +
    stat_pointinterval(aes(color=factor), scale=.75, side='left', slab_alpha=.33, show.legend=FALSE,
                       position=position_dodge(-.25, preserve='single')) +
    scale_x_discrete(name='Condition', labels=c('Reversal', 'Adversarial')) +
    ylab('Predicted Probability') +
    scale_fill_discrete(name='', labels=c('Productive\nFactor', 'Double\nPreventer',
                                          'Possible\nPreventer')) +
    facet_grid(~ vignette) +
    ggtitle('CES Sampling Probabilities') +
    ylim(0, 1) +
    theme_classic()
ggsave('plots/experiment2_ces_prob_vignette3.png', width=12, height=6)




## compare models using LOO-IC
loo_compare(list(ns=fit.ns$loo(), ces=fit.ces$loo()))

## Compare models using R-squared
R2.ns <- bayes_R2(fit.ns$draws('cause_hat', format='draws_matrix'), d$cause)
R2.ces <- bayes_R2(fit.ces$draws('cause_hat', format='draws_matrix'), d$cause)

tibble(NS=R2.ns, CES=R2.ces) %>%
    pivot_longer(NS:CES, names_to='model') %>%
    ggplot(aes(x=value, y=model)) +
    xlab('R squared') + ylab('') + xlim(0, NA) +
    stat_halfeye() +
    theme_classic()
ggsave('plots/experiment2_r2_vignette2.png', width=6, height=4)

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
ggsave('plots/experiment2_r2_ratio_vignette2.png', width=6, height=4)

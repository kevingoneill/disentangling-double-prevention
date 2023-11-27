library(tidyverse)
library(viridis)
library(scico)

## Derived model predictions
## E = PF (1 - PP) + PF PP DP

## P(E | PP) - P(E | -PP)
## P(PF)P(DP) - P(PF)

## NS = P(PP)P(E | PP) + P(-PP)P(-E | -PP, PF, DP)
##    = P(PP)P(PF)P(DP) 

d.pred <- expand_grid(p.PF=seq(0, 1, .01),
                      p.PP=seq(0, 1, .25),
                      p.DP=seq(0, 1, .01)) %>%
    group_by(p.PF, p.PP, p.DP) %>%
    mutate(p.P=p.PP * (1-p.DP),
           p.E=p.PF * (1-p.P),
           NS.PF=(1-p.PF) + p.PF*((1-p.PP)*(1-p.DP) + (1-p.PP)*p.DP + p.PP*p.DP),
           NS.DP=(1-p.DP) + p.DP*p.PF,
           NS.PP=p.PP*p.PF*p.DP,
           NS.DP_effect=NS.PF - NS.DP,
           CES.PF=sqrt((p.PF*(1-p.PF))/(p.E*(1-p.E))) * (1-p.PP + p.PP*p.DP),
           CES.DP=sqrt((p.DP*(1-p.DP))/(p.E*(1-p.E))) * p.PF * p.PP,
           CES.PP=sqrt((p.PP*(1-p.PP))/(p.E*(1-p.E))) * (p.PF*p.DP - p.PF),
           CES.DP_effect=CES.PF - CES.DP) %>%
    pivot_longer(NS.PF:CES.DP_effect, names_sep='\\.', names_to=c('model', 'variable')) %>%
    mutate(model=factor(model, levels=c('NS', 'CES')))



## Judgments of Double Preventer
d.pred %>%
    filter(variable=='DP') %>%
    ggplot(aes(x=p.PF, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_viridis(name='Causal\nJudgment of\nDouble\nPreventer', option='magma') +
    scale_x_continuous(name='P(Productive Factor)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(Double Preventer)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Double Preventer') +
    theme_classic()
ggsave('plots/dp_mean.png', width=8, height=4)




## Judgments of Productive Factor
d.pred %>%
    filter(variable=='PF') %>%
    ggplot(aes(x=p.PF, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_viridis(name='Causal\nJudgment of\nProductive\nFactor', option='magma') +
    scale_x_continuous(name='P(Productive Factor)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(Double Preventer)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Productive Factor') +
    theme_classic()
ggsave('plots/pf_mean.png', width=8, height=4)



## DP Effect (Judgments of PF - DP)
d.pred %>%
    filter(variable=='DP_effect') %>%
    ggplot(aes(x=p.PF, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_scico(name='DP Effect', palette='roma', direction=-1,
                     midpoint=0, labels=c('-1 (DP > PF)', '-.5', '0', '.5', '1 (PF > DP)')) +
    scale_x_continuous(name='P(Productive Factor)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(Double Preventer)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Double Prevention Effect') +
    theme_classic()
ggsave('plots/dp_effect.png', width=8, height=4)



## Judgments of Possible Preventer
d.pred %>%
    filter(variable=='PP') %>%
    ggplot(aes(x=p.PF, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_scico(name='Causal\nJudgment of\nPossible\nPreventer', palette='roma', direction=-1,
                     midpoint=0, labels=c('-1', '-.5', '0', '.5', '1')) +
    scale_x_continuous(name='P(Productive Factor)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(Double Preventer)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Possible Preventer') +
    theme_classic()
ggsave('plots/pp_mean.png', width=8, height=4)


d.pred %>%
    filter(p.PP %in% c(0, 1), p.DP==.1, p.PF==.99, variable %in% c('PF', 'DP')) %>%
    ggplot(aes(x=as.factor(1-p.PP), y=value, fill=variable)) +
    geom_col(position='dodge') +
    scale_x_discrete('Normality of Possible Preventer', labels=c('Normal', 'Abnormal')) +
    scale_y_continuous('Predicted Causal Judgment', breaks=c(0, .25, .5, .75, 1),
                       labels=c('0', '.25', '.5', '.75', '1')) +
    scale_fill_manual(name='', breaks=c('PF', 'DP'),
                      labels=c('Productive\nFactor\n(Normal)', 'Double\nPreventer\n(Abnormal)'),
                      values=c("#E61C38", "#C0C0C0")) +
    facet_wrap(~ model, labeller=as_labeller(c(NS='Necessity-Sufficiency',
                                               CES='Counterfactual Effect Size'))) +
    theme_classic(base_size=18)
ggsave('plots/predictions.png', width=8, height=4)

d.pred %>%
    filter(p.PP %in% c(0, 1), p.DP==.1, p.PF==.99, variable %in% c('PF', 'DP')) %>%
    ggplot(aes(x=as.factor(1-p.PP), y=value, fill=variable)) +
    geom_col(position='dodge') +
    scale_x_discrete('Normality of Possible Preventer', labels=c('Reversal', 'Adversarial')) +
    scale_y_continuous('Predicted Causal Judgment', breaks=c(0, .25, .5, .75, 1),
                       labels=c('0', '.25', '.5', '.75', '1')) +
    scale_fill_manual(name='', breaks=c('PF', 'DP'),
                      labels=c('Productive\nFactor\n(Normal)', 'Double\nPreventer\n(Abnormal)'),
                      values=c("#E61C38", "#C0C0C0")) +
    facet_wrap(~ model, labeller=as_labeller(c(NS='Necessity-Sufficiency',
                                               CES='Counterfactual Effect Size'))) +
    theme_classic(base_size=18)
ggsave('plots/predictions2.png', width=8, height=4)

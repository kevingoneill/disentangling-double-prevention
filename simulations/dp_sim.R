library(tidyverse)
library(viridis)
library(scico)

## Derived model predictions
## E = PF (1 - PP) + PF PP DP

## P(E | PP) - P(E | -PP)
## P(PF)P(DP) - P(PF)

## NS = P(PP)P(E | PP) + P(-PP)P(-E | -PP, PC, DP)
##    = P(PP)P(PF)P(DP) 

d.pred <- expand_grid(p.PC=seq(0, 1, .01),
                      p.PP=seq(0, 1, .25),
                      p.DP=seq(0, 1, .01)) %>%
    group_by(p.PC, p.PP, p.DP) %>%
    mutate(p.P=p.PP * (1-p.DP),
           p.E=p.PC * (1-p.P),
           Icard.PC=(1-p.PC) + p.PC*((1-p.PP)*(1-p.DP) + (1-p.PP)*p.DP + p.PP*p.DP),
           Icard.DP=(1-p.DP) + p.DP*p.PC,
           Icard.PP=p.PP*p.PC*p.DP,
           Icard.DP_effect=Icard.PC - Icard.DP,
           Quillien.PC=sqrt((p.PC*(1-p.PC))/(p.E*(1-p.E))) * (1-p.PP + p.PP*p.DP),
           Quillien.DP=sqrt((p.DP*(1-p.DP))/(p.E*(1-p.E))) * p.PC * p.PP,
           Quillien.PP=sqrt((p.PP*(1-p.PP))/(p.E*(1-p.E))) * (p.PC*p.DP - p.PC),
           Quillien.DP_effect=Quillien.PC - Quillien.DP) %>%
    pivot_longer(Icard.PC:Quillien.DP_effect, names_sep='\\.', names_to=c('model', 'variable'))

## Judgments of Double Preventer
d.pred %>%
    filter(variable=='DP') %>%
    ggplot(aes(x=p.PC, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_viridis(name='Causal\nJudgment of\nDouble\nPreventer', option='magma') +
    scale_x_continuous(name='P(PC)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(DP)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Double Preventer') +
    theme_classic()
ggsave('dp_mean.png', width=8, height=4)




## Judgments of Productive Cause
d.pred %>%
    filter(variable=='PC') %>%
    ggplot(aes(x=p.PC, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_viridis(name='Causal\nJudgment of\nProductive\nCause', option='magma') +
    scale_x_continuous(name='P(PC)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(DP)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Productive Cause') +
    theme_classic()
ggsave('pc_mean.png', width=8, height=4)



## DP Effect (Judgments of PC - DP)
d.pred %>%
    filter(variable=='DP_effect') %>%
    ggplot(aes(x=p.PC, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_scico(name='DP Effect', palette='roma', direction=-1,
                     midpoint=0, labels=c('-1 (DP > PC)', '-.5', '0', '.5', '1 (PC > DP)')) +
    scale_x_continuous(name='P(PC)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(DP)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Double Prevention Effect') +
    theme_classic()
ggsave('dp_effect.png', width=8, height=4)



## Judgments of Possible Preventer
d.pred %>%
    filter(variable=='PP') %>%
    ggplot(aes(x=p.PC, y=p.DP, fill=value)) +
    geom_raster() +
    facet_grid(model~p.PP, labeller=labeller(.cols=as_labeller(function(p) paste0('P(PP) = ', p)))) +
    scale_fill_scico(name='Causal\nJudgment of\nPossible\nPreventer', palette='roma', direction=-1,
                     midpoint=0, labels=c('-1', '-.5', '0', '.5', '1')) +
    scale_x_continuous(name='P(PC)', labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='P(DP)', labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(expand=FALSE) +
    ggtitle('Model Predictions: Possible Preventer') +
    theme_classic()
ggsave('pp_mean.png', width=8, height=4)






## Simulated model predictions to verify derivations
##  (this may take a while to run)
N <- 10000  ## number of samples
d.sim <- expand_grid(p.PC=seq(0, 1, .01),
            p.PP=seq(0, 1, .1),
            p.DP=seq(0, 1, .01),
            n=1:N) %>%
    group_by(p.PC, p.PP, p.DP) %>%
    mutate(PC=rbernoulli(n(), p.PC),  ## productive cause
           PP=rbernoulli(n(), p.PP),  ## possible preventer
           DP=rbernoulli(n(), p.DP),  ## double preventer
           P=PP * (1-DP),             ## actual prevention
           E=PC * (1-P)) %>%          ## actual outcome
    summarize(N=max(n),
              Quillien.PC=cor(PC, E),
              Quillien.DP=cor(DP, E),
              Icard.PC=mean((1-PC) + PC*((1-PP)*(1-DP) + (1-PP)*DP + PP*DP)),
              Icard.DP=mean((1-DP) + DP*PC))

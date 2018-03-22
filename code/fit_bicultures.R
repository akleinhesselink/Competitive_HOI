rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

fit_mod1(sim_data_file = 'data/mechanistic_sim_bicultures.rds', 
         fit_file = 'data/biculture_fit1.rds', 
         fitted_pars_file = 'data/biculture_fit1_pars.rds')

fit_mod2(sim_data_file = 'data/mechanistic_sim_bicultures.rds', 
         fit_file = 'data/biculture_fit2.rds', 
         fitted_pars_file = 'data/biculture_fit2_pars.rds')

all_fits1 <- readRDS('data/biculture_fit1.rds')
all_fits2 <- readRDS('data/biculture_fit2.rds')


plot_all_fits(all_fits1)
plot_all_fits(all_fits2)

all_fits <- left_join( all_fits1 %>% 
                         spread( model , pred_fecundity), 
                       all_fits2 %>% 
                         spread( model , pred_fecundity))
all_fits<- all_fits %>% 
  gather( model, pred_fecundity, starts_with('mod_bh') )

all_fits %>% 
  ggplot(aes( y = log(fecundity), x = log(pred_fecundity), group = model, color = focal)) + 
  geom_point() + 
  geom_smooth(aes(group = interaction(model,focal), color = focal, linetype = model), method = 'lm', se = F)




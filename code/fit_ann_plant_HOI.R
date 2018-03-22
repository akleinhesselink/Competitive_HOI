rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

fit_mod1(sim_data_file = 'data/ann_plant_sim1.rds', 
         fit_file = 'data/ann_plant_sim1_fit1.rds', 
         fitted_pars_file = 'data/ann_plant_pars_sim1_fit1.rds')

fit_mod1(sim_data_file = 'data/ann_plant_sim2.rds', 
         fit_file = 'data/ann_plant_sim2_fit1.rds', 
         fitted_pars_file = 'data/ann_plant_pars_sim2_fit1.rds')

compare_parameters( original_pars_file = 'data/ann_plant_pars1.rds', 
                    fitted_pars_file = 'data/ann_plant_pars_sim1_fit1.rds')

compare_parameters( original_pars_file = 'data/ann_plant_pars2.rds', 
                    fitted_pars_file = 'data/ann_plant_pars_sim2_fit1.rds')

all_fits <- readRDS('data/ann_plant_sim1_fit1.rds')

plot_all_fits(all_fits)



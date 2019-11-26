rm(list = ls())
source('code/plotting_parameters.R')

load('output/trade_off_results.rda')
load('output/trade_off_trials.rda')

trade_off_trials %>% 
  arrange( desc(trial), species  ) %>% 
  select(trial, species, d, m) %>% 
  mutate( round( d, 3), round(m, 3)) %>% 
  group_by( trial) %>% 
  mutate( sd(d)) 
  

plotting_curve <- 
  expand.grid( trial = 1:5, m = seq(0.05, 0.4, by = 0.01)) %>% 
  mutate( d = trade_off(m, curve_parms)) 

label_df <- data.frame( trial = 1:5, label = paste0( LETTERS[1:5], ')')) %>% 
  mutate( x_pos = -Inf, y_pos = Inf)

gg_scenarios <- 
  plotting_curve %>% 
  ggplot(aes( x = d, y = m )) + 
  geom_line() + 
  geom_point( data= trade_off_trials, aes( color = species), size = 4) + 
  geom_text( data = label_df, aes( x = x_pos, y = y_pos, label = label ), hjust = -0.5, vjust = 1.5, size = 5) + 
  facet_wrap( ~ trial, nrow = 1 ) + 
  scale_color_manual(values = my_colors[1:3], name = 'Species') + 
  ylim( c(0, 0.45)) + 
  ggtitle("Simulation Scenario") + 
  xlab( expression( 'Root tissue density, '*italic(d)*' ('*g~cm^-3*')')) +
  ylab( expression( 'Tissue respiration rate, '*italic(delta)*' ('*g~g^-1~d^-1*')')) +
  journal_theme + 
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5), 
        axis.title.y = element_text(hjust = 0.5)) + 
  theme(panel.spacing = unit(2, 'line'), 
      plot.margin = margin(1, 1, 1, 2, 'line'), 
      axis.title.x = element_text(margin = margin(1, 1, 1, 1, unit = 'line')))


ggsave(gg_scenarios , filename = 'figures/scenario_illustration_figA1.png', height = 5, width = 8, units = 'in')

tradeoff_results %>%
  filter( par_type == 'beta') %>% 
  ggplot( aes( x = trial, y = value )) + 
  geom_line() + 
  geom_point()  +
  facet_wrap(par ~ species, scales = 'free') 

tradeoff_results %>% 
  filter( par_type == 'alpha') %>%
  ggplot( aes( x = trial, y = value )) + 
  geom_line() + 
  geom_point()  +
  facet_wrap(par ~ species, scales = 'free') 



trd_off_df <- 
  tradeoff_results %>% 
  select( - par, -value ) %>% 
  distinct() %>% 
  spread( par_type, par_avg ) %>%
  mutate( HOI_ratio = beta/alpha ) %>% 
  left_join( trade_off_trials %>% group_by( trial) %>% summarise( vd = sd(d)), by = 'trial')


label_df <- 
  trd_off_df %>% 
  ungroup() %>% 
  distinct(species) %>% 
  mutate(label = paste0(LETTERS[1:3], ')')) %>% 
  mutate( x_pos = -Inf, y_pos = Inf)

gg_trd_off <- 
  trd_off_df %>% 
  ggplot( aes( x = vd, y = HOI_ratio )) + 
  geom_line() + 
  geom_point() +
  geom_text( data = label_df, aes( x = x_pos, y = y_pos , label = label), vjust = 1.5, hjust = -0.5, size = 5) + 
  facet_grid( ~ species) + 
  ylab(expression("HOI Ratio"~"("*beta*":"*alpha*")")) + 
  xlab(expression('Community Trait Difference (std. dev.'~italic(d)*')')) +
  journal_theme + 
  theme(panel.spacing = unit(2, 'line'), 
        plot.margin = margin(1, 1, 1, 2, 'line'), 
        axis.title.x = element_text(margin = margin(1, 1, 1, 1, unit = 'line')))


gg_trd_off %>% 
  ggsave(filename = 'figures/trade_off_results_figA2.png', width = 8, height = 4, units = 'in')


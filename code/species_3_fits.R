rm(list = ls())

source('code/sim_functions.R')
source('code/figure_pars.R')

nspp <- 3
start_sd <- 2
min_sd <- 1
max_refit <- 6 
data_file <- 'data/mechanistic_sim_bicultures.rds'
model <- "mod_bh2_ll"

fits <- plot <- list()

i <- 3

focal <- paste0('F', i)
  
dat <- prep_data(i, dat = readRDS(data_file))

dat12 <- 
  dat %>% 
  filter (N3 == 0)

dat13 <- 
  dat %>% 
  filter( N2 == 0)

dat23 <- 
  dat %>% 
  filter( N1 == 0)


dat_list <- list(dat12, dat13, dat23)

fit1 <- 
  lapply(dat_list, 
        FUN = fit_model, 
        form = form1,
        mod_name = model,
        start_sd = start_sd,
        min_sd = min_sd,
        max_refit = max_refit
        )


fitHOI <- 
  lapply( dat_list, 
        FUN = fit_model, 
        form = formHOI,
        mod_name = model,
        start_sd = start_sd,
        min_sd = min_sd,
        max_refit = max_refit
        )

basic_preds <- 
  mapply(x = fit1, 
       y = dat_list, 
       FUN = function(x, y, ...) predict_fit(x$par, y, ...), 
       MoreArgs = list(
         mod_name = model, 
         form = form1)
       )

HOI_preds <- 
  mapply(x = fitHOI, 
         y = dat_list, 
         FUN = function(x, y, ...) predict_fit(x$par, y, ...), 
         MoreArgs = list(
           mod_name = model, 
           form = formHOI)
         )



dat_list <- 
  mapply( x = dat_list, 
        y = as.list( as.data.frame(basic_preds) ), 
        function(x, y){ 
          x$basic <- y
          return(x) } , 
        SIMPLIFY = F
        )

dat_list <- 
  mapply( x = dat_list, 
          y = as.list( as.data.frame(HOI_preds) ), 
          function(x, y){ 
            x$HOI <- y
            return(x) } , 
          SIMPLIFY = F
  )



dat_list <- lapply( dat_list, 
        function(x) x %>% gather( type, pred, c(basic, HOI)) )

  
p1 <- 
  plot_two_sp(dat_list[[1]], focal = focal, C1 = 'N1', C2 = 'N2', C3 = 'N3') +
  geom_line(aes(y = pred, linetype = type)) +
  scale_linetype_manual(values = c(2,3)) +
  theme(legend.position = c(0.85, 0.7))
  

p2 <- 
  plot_two_sp(dat_list[[2]], focal = focal, C1 = 'N3', C2 = 'N1', C3 = 'N2') +
    geom_line(aes(y = pred, linetype = type)) +
    scale_linetype_manual(values = c(2,3)) +
    theme(legend.position = c(0.85, 0.8))
  

p3 <- 
  plot_two_sp(dat_list[[3]], focal = focal, C1 = 'N3', C2 = 'N2', C3 = 'N1') +
    geom_line(aes(y = pred, linetype = type)) +
    scale_linetype_manual(values = c(2,3)) +
    theme(legend.position = c(0.7, 0.75), legend.box.just = 'right')

p1
p2
p3  





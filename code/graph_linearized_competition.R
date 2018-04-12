
data_file <- 'data/mechanistic_sim_bicultures.rds'

dat <- readRDS(data_file)

test <- 
  dat %>% 
  group_by( focal )  %>% 
  mutate( f2 = max(fecundity)/fecundity ) %>% 
  filter( focal == 'F1') %>% 
  filter( N1 == 0)

test %>% 
  ggplot(aes( x = N2, y = N3, z = f2)) + 
  geom_contour()


test <- 
  dat %>% 
  group_by( focal )  %>% 
  mutate( f2 = max(fecundity)/fecundity ) %>% 
  filter( focal == 'F2') %>% 
  filter( N2 == 0)

test %>% 
  ggplot(aes( x = N1, y = N3, z = f2)) + 
  geom_contour()


test <- 
  dat %>% 
  group_by( focal )  %>% 
  mutate( f2 = max(fecundity)/fecundity ) %>% 
  filter( focal == 'F3') %>% 
  filter( N2 == 0)

test %>% 
  ggplot(aes( x = N3, y = N1, z = f2)) + 
  geom_contour()

m1 <- lm(f2 ~ N1 + N3 + I(N1*N3), data = test)
summary(m1)

gam1 <- mgcv::gam(f2 ~ N1 + N3 + s(I(N1*N3)), data = test)
summary(m1)

test$m1 <- predict(m1)
test$gam1 <- predict(gam1)

test %>% 
  ggplot( aes( x = N3, y = N1, z = f2)) + 
  geom_contour()

test %>% 
  ggplot( aes( x = N3, y = N1, z = m1)) + 
  geom_contour()

test %>% 
  ggplot( aes( x = N3, y = N1, z = gam1)) + 
  geom_point() + 
  geom_contour()

plot(gam1)


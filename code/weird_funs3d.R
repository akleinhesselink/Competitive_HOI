rm(list = ls())
library(plotly)

Sys.setenv("plotly_username"="arklein")
Sys.setenv("plotly_api_key"= "zLHQXaXsi5Bpx4ZnFOwD")

library(tidyverse)

weird_fun <- function( x, parms ) { 
  
  with(parms, { 
    
    r*(1 - a*x[, 1] )*(1 - b*x[, 2])  
  })
}


mypars <- list( a = 0.3, b = 0.1, r = 0.2 )

x <- expand.grid( seq( 0,  10, by = 0.1), seq(0, 10, by = 0.1) )

x$z <- weird_fun(x, mypars)

test <- data.frame (x) %>% 
 spread( Var2, z) %>% as.matrix 

weird1 <- test[, -1]

weird_fun2 <- function( x, parms) { 
  
  with(parms, { 
    r*(1/(1 + a*x[, 1]))*(1/(1 + b*x[, 2]))  
  })
}
  
  
x$z <-  weird_fun2(x, mypars )

test <- data.frame (x) %>% 
  spread( Var2, z) %>% as.matrix 

weird2 <- test[, -1]

make_3d_plot <- function(surface, my_title) { 
  
  p <- plot_ly(z = ~surface ) %>% 
    add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      )
    ) %>%
    layout(
      title = my_title,
      scene = list(
        camera=list(
          eye = list(x=1.87, y=0.88, z=0.5)
        )
      )
    )
  return( p )
}

p1 <- make_3d_plot(test, 'weird1')

link1 <- api_create(p1, filename = "weird1")

p2 <- make_3d_plot(weird2, 'weird2')
link2 <- api_create(p2, filename = "weird2")


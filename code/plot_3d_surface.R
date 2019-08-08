rm(list = ls())
library(plotly)

Sys.setenv("plotly_username"="arklein")
Sys.setenv("plotly_api_key"= "zLHQXaXsi5Bpx4ZnFOwD")

load('output/pred_surfaces.rda')

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

p1 <- make_3d_plot(pred_surface1, 'species 1 early')
p2 <- make_3d_plot(pred_surface2, 'species 2 middle')
p3 <- make_3d_plot(pred_surface3, 'species 3 late')

p4 <- make_3d_plot( pred_surface4, 'species 1 test')


link1 <- api_create(p1, filename = "species1_surface")
link2 <- api_create(p2, filename = "species2_surface")
link3 <- api_create(p3, filename = "species3_surface")



